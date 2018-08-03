import config
import csv
import datetime as dt
import re
from point import Point
from episode import Episode
from location import Location
from misc_funcs import (distance, inner_angle_sphere, kde, min_peak, 
	state_transition_matrix, viterbi, emission_probabilities)
from math import sqrt

class Trace(object):
	"""
	A "trace", a GPS trace, is all the data associated with one itinerum user.
	It's mainly treated here as a temporal/spatial sequence of points.
	"""

	def __init__(self, user_id, raw_data, raw_survey):
		"""
		Construct the user object by pulling all data pertaining to this user, 
		identified by ID.
		"""
		self.id = user_id
		self.raw = raw_data
		self.home = raw_survey[0]  # Raw survey data passed as a list of 3 locations
		self.work = raw_survey[1]
		self.school = raw_survey[2]
		# the set of original input points, minus any that get cleaned out. 
		self.points = []
		# points removed during cleaning
		self.discarded_points = []
		# "known_subsets" is a list of lists of points partitioned by gaps in
		# the data where e.g. the phone has turned off inexplicably.
		self.known_subsets = []
		# "known_subsets_interpolated" are the subsets with points added in the
		# middle. "all_interpolated_points" is the flattened version of the
		# preceding, containing all real and interpolated points in one place.
		# This one gets used for KDE etc.
		self.known_subsets_interpolated = []
		# list of potential activity locations
		self.locations = []
		# dictionary of activity episodes?
		self.episodes = []
		# records the number of activity records written so far
		self.activity_count = 0
		# number of identical timestamps for this user
		self.identical = 0
		# read in all time and location data for points
		# right now only using a few fields
		for row in raw_data:
			self.points.append( Point(
				int(row['timestamp_epoch']),
				float(row['longitude']), float(row['latitude']), 
				float(row['h_accuracy'])
			))
		# sort the list by time
		self.points.sort(key=lambda x: x.unix_time)
		# measure to and from neighbors
		all_indices = [i for i, p in enumerate(self.points)]
		self.observe_neighbors(all_indices)

	@property
	def all_interpolated_points(self):
		"""Get all (real & interpolated) points in one big list"""
		return [ p for s in self.known_subsets_interpolated for p in s ]

	def interpolate_segment(self, segment):
		"""
		Takes a known subset (a list of ordered points) and interpolates
		linearly between them such that the distance between new points
		is never greater than a value specified in config, e.g. 30m.
		"""
		new_points = []
		# For each point bu the last
		for i in range(0, len(segment) - 1):
			pair_int = segment[i].pair_interpolation(segment[i+1])
			new_points.extend(pair_int)
		new_points.append(segment[-1])
		return new_points

	def make_known_subsets(self):
		"""
		Partition the trace points into contiguous sets for which we're confident
		we don't have substantial missing data. That is, exclude segments where it 
		seems like we have no data, but substantial movement; for which trip and 
		activity reconstruction would be impossible.
		TODO: Eventually we will need to check for subway trips.
		"""
		segment = [self.points[0]]
		# iterate over all points (except the first). Test each point to see
		# whether we add it to the current segment or the one after.
		for i in range(1, len(self.points)):
			if distance(self.points[i], self.points[i-1]) > 1000:
				# append point to next segment
				self.known_subsets.append(segment)
				segment = [self.points[i]]
			else:
				# add this point to the current segment
				segment.append(self.points[i])
		self.known_subsets.append(segment)
		# check these segments for sufficient length
		self.known_subsets = [ s for s in self.known_subsets if len(s) >=5 ]
		self.known_subsets = [ s for s in self.known_subsets if s[-1].unix_time - s[0].unix_time > 3600 ]
		for seg_index, segment in enumerate(self.known_subsets):
			for point in segment:
				point.known_subset = seg_index
		if config.db_out:
			print('\t', len(self.known_subsets) - 1, 'gap(s) found in data')

	def get_activity_locations(self):
		"""
		Get activity locations for this trace. (Create inputs for a KDE
		function and find peaks in the surface.)
		"""
		for subset in self.known_subsets:
			# interpolate the subset and weight the points
			interpolated_subset = subset  # self.interpolate_segment(subset)
			self.known_subsets_interpolated.append(interpolated_subset)
			self.weight_points(interpolated_subset)
		if len(self.all_interpolated_points) > 75000:
			raise Exception('Too many points for efficient KDE')
		# format as vectors for KDE function
		Xvector = [p.x for p in self.all_interpolated_points]
		Yvector = [p.y for p in self.all_interpolated_points]
		Wvector = [p.weight for p in self.all_interpolated_points]
		# run the KDE, returning density estimates at the input points
		estimates = kde(Xvector, Yvector, Wvector)
		# assign probability estimates to points
		for point, prob in zip(self.all_interpolated_points, estimates):
			point.kde_p = prob
		# determine average GPS accuracy value for this user
		# (sqrt of the mean variance)
		mean_accuracy = sqrt(sum([p.accuracy**2 for p in self.points]) / len(self.points))
		# estimate peak threshold value
		threshold = min_peak(
			mean_accuracy,  # mean sd of GPS accuracy for user
			sum(Wvector)  # total seconds entering KDE
		) 
		# Find peaks in the density surface
		locations = self.find_peaks(threshold)
		# store the result
		locations = self.locations.extend(locations)
		return self.locations

	def identify_locations(self):
		"""
		Identify locations with user-provided home, work, school locations if
		possible. This algorithm was written in a hurry and needs to be made
		much more robust. It is not currently called anywhere in the code. TODO
		"""
		if self.home:
			for location in self.locations:
				if distance(location, self.home) <= 150:  # meters
					location.identify('home')
		if self.work:
			for location in self.locations:
				if distance(location, self.work) <= 150:  # meters
					location.identify('work')
		if self.school:
			for location in self.locations:
				if distance(location, self.school) <= 150:  # meters
					location.identify('school')

	def break_trips(self):
		"""
		Use a Hidden Markov Model to classify all points as deriving from
		either time spent travelling or time spent at one of the potential
		activity locations. Allocate time to these sequences of activities
		accordingly.
		"""
		# get list of potential state indices
		# 0 is travel, others are then +1 from their list location
		states = range(0, len(self.locations) + 1)
		# initial_state probabilities: 50/50 travelling or stationary
		start_probs = [0.5] + [(0.5/len(states))] * len(states)
		# list of locations that actually get used
		used_locations = set()
		# run the viterbi algorithm on each known subset
		for points in self.known_subsets_interpolated:
			emission_probs = []
			for point in points:
				emission_probs.append(emission_probabilities(point, self.locations))
			# run viterbi on each subset
			state_path = viterbi(states, emission_probs, start_probs,
					state_transition_matrix(states))
			# store state classification and location references with points
			for i, state in enumerate(state_path):
				points[i].state = state
				points[i].location = None if state == 0 else self.locations[state-1]
			# note which locations have actually been used
			for location in set([p.location for p in points if p.location]):
				location.visited = True
			# iterate over points/states, creating episodes at transitions
			for i, point in enumerate(points):
				if i == 0:
					# record first episode
					self.episodes.append( Episode(point.time, point.location) )
				else:
					# look for state changes
					if prev_point.state != point.state:
						# assume the transition happens halfway between points
						transition_time = prev_point.time+(point.time-prev_point.time)/2
						self.episodes.append( Episode(transition_time, point.location) )
				prev_point = point
			# unknown time ends every known segment
			self.episodes.append( Episode(points[-1].time, is_unknown_time=True) )
		if config.db_out:
			print('\tFound', len(self.episodes), 'episodes')

	def find_peaks(self, threshold):
		"""
		Detect peaks in the KDE surface which are above the time-spent
		threshold. KDE values are stored in self.points.
		"""
		points = [point for point in self.all_interpolated_points if point.kde_p >= threshold]
		if config.db_out:
			print('\tClustering', len(points), 'points above', threshold, 'threshold')
		# For each point:
		#   for every other point within cluster distance:
		#      if comparison point has higher KDE value, this is not the peak
		potential_activity_locations = []
		loc_num = 1
		for point in points:
			is_peak = True  # starting assumption
			for neighbor in points:
				if distance(point, neighbor) > config.location_distance:
					continue
				if point.kde_p < neighbor.kde_p:
					is_peak = False  # assumption proven false if anything else higher
					break
			if is_peak and point not in potential_activity_locations:
				location = Location(point.longitude, point.latitude, loc_num)
				potential_activity_locations.append(location)
				loc_num += 1
		return potential_activity_locations

	def weight_points(self,segment):
		"""
		Assign time-based weights to a series of sequential points.
		Values are in seconds, and are split between neighboring points.
		e.g.  |--p2-time---|
		p1----|----p2------|------p3

		"""
		assert len(segment) > 1
		# iterate over all points except terminal.
		for i in range(len(segment)-1):
			this_point = segment[i]
			next_point = segment[i+1]

			d = distance(this_point, next_point)
			t = (next_point.time - this_point.time).total_seconds()
			if t == 0:  # occurs when there are identical timestamps
				self.identical += 1  # for accounting 
				# there is no time weight between these points
				w1 = 0
				w2 = 0
			else:
				w1 = this_point.weight_decimal( config.weight_coef * d / t ) * t
				w2 = ( 1 - this_point.weight_decimal( config.weight_coef * d / t) ) * t
			this_point.add_weight(w1)
			next_point.add_weight(w2)

	def remove_repeated_points(self):
		"""There are some records in the coordinates table that are simply
		reapeted verbatim. Points are already sorted by time, so to find
		these we just need to loop through and compare adjacent points.
		"""
		unique_points = []
		to_remove = []
		for i, point in enumerate(self.points):
			# TODO can we do this better with an __eq__ function?
			uid = str(point.latitude)+str(point.longitude)+str(point.time)
			# if this is the first we've seen this exact record
			if uid not in unique_points:
				unique_points.append(uid)
			else: # we've already seen this exact point
				to_remove.append(i)
		# remove the points from the main list to the recycling bin
		for i in reversed(to_remove):
			self.pop_point(i)
		if config.db_out and len(to_remove) > 0:
			print('\t', len(to_remove), 'points removed as exact duplicate')

	def pop_point(self, key):
		"""
		Pop a point off the current list and add it to the discard bin.
		Then update it's former neighbors in the list.
		"""
		# pop and append
		point = self.points.pop(key)
		point.discarded = True
		self.discarded_points.append(point)
		# now using the (former) key, update the (former) neighbors
		i_ante = key - 1
		i_post = key  # has already shifted over
		self.observe_neighbors([i_ante, i_post])

	def find_duplicates(self):
		"""
		Find any repeated point locations.
		"""
		# dictionary of unique locations, with a coordinate string as key
		# containing a list of points with that exact location
		locations = {}
		for point in self.points:
			key = str(point.geom)
			if key not in locations:
				locations[key] = [point]
			else:
				locations[key].append(point)
		for k in locations.keys():
			if len(locations[k]) > 1:
				if config.db_out:
					print(k, locations[k])

	def observe_neighbors(self, indices=[]):
		"""
		Get angle and distance measurements to and from adjacent points.
		Store in the point object. Operates on points, the inices of which
		are given in a list referring to the current point list.
		"""
		# for each point, except those first and last
		for i in indices:
			# skip first and last points
			if i <= 0 or i >= len(self.points)-1:
				continue
			point = self.points[i]
			# Find the nearest DIFFERENT geometries
			# first walk backwards to the next different point
			i_ante = i-1
			while point.geom == self.points[i_ante].geom:
				if i_ante <= 0:
					break
				i_ante -= 1
			# then forward to the next different point
			i_post = i+1
			while point.geom == self.points[i_post].geom:
				if i_post >= len(self.points) - 1:
					break
				i_post += 1
			# distance to previous point
			point.d_ante = distance(point, self.points[i_ante])
			# distance to next point
			point.d_post = distance(point, self.points[i_post])
			# if either distance is zero, this is an end point and has no angle
			if point.d_ante * point.d_post == 0:
				point.angle = 180
			else:
				# calculate the inner angle
				point.angle = inner_angle_sphere(self.points[i_ante],
					point, self.points[i_post])
			# is point is identical with both neighbors?
			if (point.geom == self.points[i-1].geom and
				point.geom == self.points[i+1].geom):
				point.inter = True

	def remove_sequential_duplicates(self):
		"""
		Remove points where both neighbors have identical locations.
		"""
		to_remove = []
		for i, point in enumerate(self.points):
			# if there is no distance between this and neighboring points
			if point.inter:
				to_remove.append(i)
		# remove the points from the main list to the recycling bin
		for i in reversed(to_remove):
			self.pop_point(i)
		if config.db_out and len(to_remove) > 0:
			print('\t', len(to_remove), 'points removed as duplicate')

	def remove_known_error(self, error_limit):
		"""
		Remove points reporting positional error beyond a given limit (meters).
		"""
		to_remove = []
		for i, point in enumerate(self.points):
			if point.accuracy > error_limit:
				to_remove.append(i)
		# remove the points from the main list to the recycling bin
		for i in reversed(to_remove):
			self.pop_point(i)
		if config.db_out and len(to_remove) > 0:
			print('\t', len(to_remove), 'points removed as high stated error')

	def remove_positional_error(self):
		"""
		Use angle/distance based cleaning rules.
		"""
		i = self.find_error_index()
		count = 0
		while i:
			self.pop_point(i)
			i = self.find_error_index()
			count += 1
		if config.db_out and count > 0:
			print('\t', count, 'points removed by positional cleaning')

	def find_error_index(self):
		"""
		Returns the index of the craziest point.
		"""
		# check first for angle == 0, as these are all obviously crazy
		for i, point in enumerate(self.points):
			if point.angle == 0:
				return i
		# check for very low angle, moderately high distance
		# and return the craziest first
		errors = {}
		for i, point in enumerate(self.points):
			# all but first and last
			if i == 0 or i == len(self.points) - 1:
				continue
			nearest = min(point.d_ante, point.d_post)
			if point.angle < 15 and nearest > 100:
				error_index = 1 / (point.angle / nearest)
				errors[error_index] = i
			if len(errors.keys()) > 0:
				return errors[max(errors.keys())]
		return False

	def get_days(self):
		"""
		Slice episodes by day.
		"""
		days = {}
		# for each episode except the last (always compare to next)
		for i in range(0, len(self.episodes) - 1):
			# what date(s) did this occur on?
			tz = self.episodes[i].start.tzinfo
			start_date = self.episodes[i].start.date()
			end_date = self.episodes[i+1].start.date()
			date_range = (end_date-start_date).days
			# for each date touched by each activity
			for offset in range(0, date_range + 1):
				date = start_date + dt.timedelta(days=offset)
				# initialize the date the first time we see it
				if date not in days:
					days[date] = {'home': [], 'work': [], 'school': [],
					'other': [], 'travel': [], 'unknown': [], 'total': []}
				# how much of this activity occured on this date?
				# first limit start_time to start of this day
				slicer1 = dt.datetime.combine(date, dt.time(0,0,0,0,tzinfo=tz) )
				st1 = slicer1 if self.episodes[i].start < slicer1 else\
					self.episodes[i].start
				# now limit end time to end of this day
				slicer2 = dt.datetime.combine(date + dt.timedelta(days=1), dt.time(0,0,0,0,tzinfo=tz) )
				st2 = slicer2 if self.episodes[i+1].start > slicer2 else\
					self.episodes[i+1].start
				# get duration in minutes from timedelta obj
				duration = (st2-st1).total_seconds() / 60
				# add activity duration to the total for this date
				days[date]['total'].append(duration)
				# now also add the duration to the appropriate category
				days[date][self.episodes[i].type].append(duration)
		return days

	def flush(self):
		"""
		After everything is finished, write all the output from this trace to CSV 
		files defined in config. All writing to files should be done here if 
		possible. Any data that needs to ultimately find it's way here should be 
		stored as a property. All Trace objects call this at the end, and the 
		files are initialized in main, so we only append rows here.
		"""
		self.write_locations()
		self.write_episodes()
		self.write_points()
		self.write_day_summary()

	def write_locations(self):
		""" Output activity locations to CSV."""
		# write potential activity locations to file
		with open(config.output_locations_file, "a") as f:
			for location in self.locations:
				f.write("{},{},{},{},{},{}\n".format(
					self.id,           # user_id
					location.id,       # location_id
					location.longitude,
					location.latitude,
					location.name,     # description
					location.visited   # whether location was used or not
				))  

	def write_episodes(self):
		""" Output episode data to CSV."""
		# write episodes file
		with open(config.output_episodes_file, "a") as f:
			for i, episode in enumerate(self.episodes):
				f.write("{},{},{},{},{},{},{}\n".format(
					self.id,  # user_id
					i,        # activity sequence
					episode.location_id if episode.location else '',
					'',       # mode (not currently used)
					True if episode.unknown else '',
					episode.start,             # timestamp
					episode.start.timestamp()  # unix time
			)) 

	def write_points(self):
		""" Output point attributes to CSV for debugging."""
		# write preliminary points file
		# 'user_id,lon,lat,removed,interpolated,state'
		with open(config.output_points_file, 'a') as f:
			for point in self.all_interpolated_points + self.discarded_points:
				f.write("{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
					self.id,
					point.unix_time,
					point.known_subset if point.known_subset is not None else '',
					point.longitude,
					point.latitude,
					point.x,
					point.y,
					point.weight,
					point.discarded,
					point.synthetic,
					point.state if point.state is not None else '',
					point.kde_p if point.kde_p is not None else ''))

	def write_day_summary(self):
		"""Output daily summary to CSV."""
		days = self.get_days()
		with open(config.output_days_file, 'a') as f:
			for date, data in days.items():
				f.write("{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
					self.id, date, date.weekday(),
					sum(data['total']),
					len(data['travel']),
					sum(data['travel']),
					sum(data['unknown']),
					sum(data['home']),
					sum(data['work']),
					sum(data['school']),
					len(data['home']),
					len(data['work']),
					len(data['school'])
				))
