import config
import csv
import datetime as dt
from points import GPSpoint, Location
from episode import Episode
from HMM import viterbi, state_transition_matrix, emission_probabilities
from gaussian import kde, min_peak
from spatial_functions import project, unproject, distance, inner_angle_sphere
from math import sqrt, ceil


class Trace(object):
	"""
	A "trace", a GPS trace, is all the data associated with one itinerum user.
	It's mainly treated here as a temporal/spatial sequence of points.
	"""

	def __init__(self, user_id, raw_data, named_places, locks):
		"""
		Construct the user object by pulling all data pertaining to this user,
		identified by ID.
		"""
		self.locks = locks
		self.id = user_id
		self.raw = raw_data
		self.named_places = named_places
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
		# list of activity episodes (travel, activity, home, work, etc)
		self.episodes = []
		# records the number of activity records written so far
		self.activity_count = 0
		# number of identical timestamps for this user
		self.identical = 0
		# read in all time and location data for points
		# right now only using a few fields
		for row in raw_data:
			self.points.append(GPSpoint(int(row['timestamp_epoch']),
			float(row['longitude']), float(row['latitude']), float(row['h_accuracy'])))
		# sort the list by time
		self.points.sort(key=lambda x: x.unix_time)
		# measure to and from neighbors
		all_indices = [i for i, p in enumerate(self.points)]
		self.observe_neighbors(all_indices)

	@property
	def all_interpolated_points(self):
		"""Get all (real & interpolated) points in one big list"""
		return [p for s in self.known_subsets_interpolated for p in s]

	def temporally_interpolate_points(self, segment):
		"""
		This function interpolates linearly in spacetime where the temporal
		gap between two points is sufficiently large.
		"""
		max_time_gap = 600  # seconds
		new_points = []
		# For each segment (pair of points)
		for i in range(1, len(segment)):
			p1, p2 = segment[i-1], segment[i]
			delta_t = p1.delta_t(p2)
			if delta_t <= max_time_gap:
				new_points.append(p1)  # no interpolation to do
				continue
			# spatial
			dist = distance(p1, p2)
			n_segs = ceil(delta_t / max_time_gap)
			seg_span = delta_t / n_segs
			delta_x = (p2.x - p1.x) / n_segs
			delta_y = (p2.y - p1.y) / n_segs
			# iteration over new subsegments
			inter_points = []
			acc = (p1.accuracy + p2.accuracy) / 2
			for j in range(1, n_segs):
				lng, lat = unproject(p1.x + j * delta_x, p1.y + j * delta_y)
				t = p1.unix_time + j * delta_t / n_segs
				new_point = GPSpoint(t, lng, lat, acc)
				new_point.synthetic = True
				inter_points.append(new_point)
			new_points.extend([p1] + inter_points)
		# add the last point
		new_points.append(segment[-1])
		return new_points

	def spatially_interpolate_points(self, segment):
		"""
		Takes a list of ordered points and interpolates spatially between them
		such that the distance between the returned list of points is never
		greater than a value specified in config, e.g. 30m.
		Temporally, there are two paradigms. For segments faster than walking
		speed, time is allocated uniformly. For those slower, it is assumed that
		walking-speed travel begins at the last possible moment, allowing time to
		accumulate at the preceding point. This function does NOT assign temporal
		weights; those are applied later according to timestamps.
		"""
		new_points = []
		# For each segment
		for i in range(1, len(segment)):
			p1, p2 = segment[i-1], segment[i]
			dist = distance(p1, p2)
			if dist <= config.interpolation_distance:
				new_points.append(p1)  # no interpolation to do
				continue
			walk_speed = 5*1000/3600  # 5 kmph in mps
			delta_t = p1.delta_t(p2)
			# spatial
			n_segs = ceil(dist / config.interpolation_distance)
			seg_len = dist / n_segs
			delta_x = (p2.x - p1.x) / n_segs
			delta_y = (p2.y - p1.y) / n_segs
			# iteration over new subsegments
			inter_points = []
			acc = (p1.accuracy + p2.accuracy) / 2
			for j in range(1, n_segs):
				lng, lat = unproject(p1.x + j * delta_x, p1.y + j * delta_y)
				if dist >= delta_t * walk_speed:
					# if faster than walking speed, assign time uniformly
					t = p1.unix_time + j * delta_t / n_segs
				else:
					# slower than walking speed, so assign time backwards from last
					# point at walking speed
					t = p2.unix_time - (dist/n_segs)/walk_speed*(n_segs-j)
				new_point = GPSpoint(t, lng, lat, acc)
				new_point.synthetic = True
				inter_points.append(new_point)
			new_points.extend([p1] + inter_points)
		# add the last point
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
		self.known_subsets = [s for s in self.known_subsets if len(s) >= 5]
		self.known_subsets = [s for s in self.known_subsets if
					s[-1].delta_t(s[0]) > 3600]
		for seg_index, segment in enumerate(self.known_subsets):
			for point in segment:
				point.known_subset = seg_index
		if config.debug_output:
			print('\t', len(self.known_subsets) - 1, 'gap(s) found in data')

#	def partial_interpolation_removal(self, segment, cutoff):
#		# TODO this function causes an error (sequential locations for user A)
#		# TODO needs documentation
#		new_segment = []
#		non_synth = [segment.index(point)
#				for point in segment if not point.synthetic]
#		for i in range(len(non_synth) - 1):
#			if non_synth[i+1] - non_synth[i] > cutoff:
#				# keep first endpoint
#				# second endpoint appended on next iteration
#				new_segment.append(segment[non_synth[i]])
#			else:
#				# keep everything but the last endpoint
#				sub_segment = [segment[j]terpolate
#					for j in range(non_synth[i], non_synth[i+1])]
#				new_segment.extend(sub_segment)
#		new_segment.append(segment[non_synth[-1]])
#		return new_segment

	def get_activity_locations(self):
		"""
		Get activity locations for this trace. (Create inputs for a KDE
		function and find peaks in the surface.)
		"""
		for subset in self.known_subsets:
			# interpolate the subset and weight the points
			interpolated_subset = self.spatially_interpolate_points(subset)
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
		mean_accuracy = sqrt(sum([p.accuracy**2 for p in self.points]) /
			len(self.points))
		# estimate peak threshold value
		threshold = min_peak(mean_accuracy,  # mean sd of GPS accuracy for user
				sum(Wvector))  # total seconds entering KDE
		# Find peaks in the density surface
		locations = self.find_peaks(threshold)
		# store the result
		locations = self.locations.extend(locations)
		return self.locations

	def identify_named_locations(self):
		"""
		Take user supplied (named) places and apply the labels to discovered
		locations if possible. This is based simply on distance, though this
		function is called after time has been allocated, so we could in theory
		base it on a function of distance and time spent.
		"""
		# for each named place, check distance to other locations
		for name, place in self.named_places.items():
			dists = [distance(place, loc) for loc in self.locations]
			if min(dists) <= 200:  # meters
				i = dists.index(min(dists))
				self.locations[i].identify(name)

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
		# do temporal interpolation on all known subsets
		self.known_subsets_interpolated = [self.temporally_interpolate_points(subset)
			for subset in self.known_subsets_interpolated]
		# run the viterbi algorithm on each known subset
		for points in self.known_subsets_interpolated:
			emission_probs = emission_probabilities(points, self.locations)
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
					self.episodes.append(Episode(point.time, point.location))
				else:
					# look for state changes
					if prev_point.state != point.state:
						# assume the transition happens halfway between points
						transition_time = prev_point.time+(point.time-prev_point.time)/2
						self.episodes.append(Episode(transition_time, point.location))
				prev_point = point
			# unknown time ends every known segment
			self.episodes.append(Episode(points[-1].time, is_unknown_time=True))
		if config.debug_output:
			print('\tFound', len(self.episodes), 'episodes')

	def find_peaks(self, threshold):
		"""
		Detect peaks in the KDE surface which are above the time-spent
		threshold. KDE values are stored in self.points.
		"""
		points = [point for point in self.all_interpolated_points
			if point.kde_p >= threshold]
		if config.debug_output:
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

	def weight_points(self, segment):
		"""
		Assign time-based weights to a series of sequential points.
		Values are in seconds, and split evenly between neighboring points, e.g.:
		 |-p1-time-|--p2-time-|...etc
		p1---------|---------p2------|------p3
				      ^midpoint
		"""
		assert len(segment) > 1
		# iterate over all sub-segments
		for i in range(len(segment)-1):
			p1, p2 = segment[i], segment[i+1]
			time_diff = p1.delta_t(p2)  # time in seconds
			p1.add_weight(time_diff/2)
			p2.add_weight(time_diff/2)

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
			else:  # we've already seen this exact point
				to_remove.append(i)
		# remove the points from the main list to the recycling bin
		for i in reversed(to_remove):
			self.pop_point(i)
		if config.debug_output and len(to_remove) > 0:
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
				if config.debug_output:
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
		if config.debug_output and len(to_remove) > 0:
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
		if config.debug_output and len(to_remove) > 0:
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
		if config.debug_output and count > 0:
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
				slicer1 = dt.datetime.combine(date, dt.time(0, 0, 0, 0, tzinfo=tz))
				st1 = slicer1 if self.episodes[i].start < slicer1 else\
					self.episodes[i].start
				# now limit end time to end of this day
				slicer2 = dt.datetime.combine(date + dt.timedelta(days=1),
					dt.time(0, 0, 0, 0, tzinfo=tz))
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
		while self.remove_short_stationary_episodes():
			pass
		self.write_locations()
		self.write_episodes()
		self.write_points()
		self.write_day_summary()

	def write_locations(self):
		""" Output activity locations to CSV."""
		# write potential activity locations to file
		with open(config.output_dir+'/locations.csv', "a") as f:
			for location in self.locations:
				if config.multi_process:
					self.locks[0].acquire()
				f.write("{},{},{},{},{},{}\n".format(self.id,  # user_id
					location.id,       # location_id
					location.longitude,
					location.latitude,
					location.name,     # description
					location.visited))   # whether location was used or not
				if config.multi_process:
					self.locks[0].release()

	def write_episodes(self):
		""" Output episode data to CSV."""
		# write episodes file
		with open(config.output_dir+'/episodes.csv', "a") as f:
			for i, episode in enumerate(self.episodes):
				attributes = [self.id, i,  # activity sequence
					episode.location_id if episode.location else '', '',
					# mode (not currently used)
					True if episode.unknown else '',
					episode.start,               # timestamp
					episode.start.timestamp()]  # unix time
				if config.multi_process:
					self.locks[1].acquire()
				f.write(','.join([str(a) for a in attributes])+'\n')
				if config.multi_process:
					self.locks[1].release()

	def write_points(self):
		""" Output point attributes to CSV for debugging."""
		# write preliminary points file
		# 'user_id,lon,lat,removed,interpolated,state'
		with open(config.output_dir+'/classified_points.csv', 'a') as f:
			for point in self.all_interpolated_points + self.discarded_points:
				attributes = [self.id, point.unix_time,
					point.known_subset if point.known_subset is not None else '',
					point.longitude, point.latitude, point.x, point.y,
					point.weight, point.discarded,
					point.human_timestamp, point.synthetic,
					point.state if point.state is not None else '',
					point.kde_p if point.kde_p is not None else '']
				if config.multi_process:
					self.locks[2].acquire()
				f.write(','.join([str(a) for a in attributes])+'\n')
				if config.multi_process:
					self.locks[2].release()

	def write_day_summary(self):
		"""Output daily summary to CSV."""
		days = self.get_days()
		with open(config.output_dir+'/days.csv', 'a') as f:
			for date, data in days.items():
				attributes = [self.id, date, date.weekday(), sum(data['total']),
len(data['travel']), sum(data['travel']), sum(data['unknown']),
sum(data['home']), sum(data['work']), sum(data['school']),
len(data['home']), len(data['work']), len(data['school'])]
				if config.multi_process:
					self.locks[3].acquire()
				f.write(','.join([str(a) for a in attributes])+'\n')
				if config.multi_process:
					self.locks[3].release()

	def remove_short_stationary_episodes(self):
		"""
		Identify and remove stationary episodes shorter than a given time.
		"""
		self.episodes.sort()
		for i in range(0, len(self.episodes)-1):
			cur, nxt = self.episodes[i], self.episodes[i+1]
			# we don't care how short travel and unknown episodes are
			if cur.type in ['unknown', 'travel']:
				continue
			delta_seconds = (nxt.start - cur.start).total_seconds()
			# if too short
			if delta_seconds < config.minimum_activity_time:
				self.remove_episode(i)
				return True
		return False

	def remove_episode(self, i):
		"""
		Given a stationary episode, identified by index position, remove the
		episode and dissolve its time into surrounding episodes as appropriate.
		"""
		assert 0 <= i <= len(self.episodes)
		prev = self.episodes[i-1] if i-1 >= 0 else None
		this = self.episodes[i]
		nxt = self.episodes[i+1] if i+1 <= len(self.episodes) else None
		assert this.type not in ['unknown', 'travel']
		if nxt and prev:
			if prev.type == nxt.type:
				# these two get dissolved into prev
				self.episodes.pop(i+1)  # pop next
				self.episodes.pop(i)   # pop this
			else:
				# draw the times from the surrounding episodes into the middle
				# end_prev --> deleted_episode <-- start_next
				# start next is the only one that actually needs updated
				nxt.start = nxt.start + (nxt.start - this.start)/2
				self.episodes.pop(i)
		else:  # is at either end of the sequence
			self.episodes.pop(i)
