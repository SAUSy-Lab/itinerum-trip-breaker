import config, csv
from point import Point
from episode import Episode
from location import Location
from misc_funcs import distance, inner_angle_sphere, kde, min_peak, gaussian, ts_str, unproject
from datetime import timedelta, datetime
from random import sample
from math import sqrt

class Trace(object):
	"""A "trace", a GPS trace, is all the data associated with one itinerum user.
		It's mainly treated here as a temporal/spatial sequence of points."""

	def __init__(self,user_id, raw_data): #TODO big ol bottleneck
		"""Construct the user object by pulling all data pertaining to this user.
			Identified by ID"""
		self.id = user_id		# 
		# There are many lists of points. The first is the set of original
		# input points, minus any that get cleaned out. Those are moved to 
		# "discarded_points" where they are left to their own devices. 
		# "known_subsets" is a list of lists of points partitioned by gaps in 
		# the data where e.g. the phone has turned off inexplicably.
		# "known_subsets_interpolated" are the subsets with points added in the 
		# middle. "all_interpolated_points" is the flattened version of the 
		# preceding, containing all real and interpolated points in one place.
		# This one gets used for KDE etc.
                self.raw = raw_data
		self.points = []
		self.discarded_points = []
		self.known_subsets = []
		self.known_subsets_interpolated = []
		self.all_interpolated_points = []
		# list of potential activity locations
		self.locations = []
		# dictionary of activity episodes?
		self.episodes = []
		# records the number of activity records written so far
		self.activity_count = 0
		# home location
		self.home = None
		# work location
		self.work = None
		# school location 
		self.school = None
		# read in all time and location data for points
		# right now only using a few fields

		for row in raw_data:
			self.points.append(
				Point( #works for our format not others 
					row[7],
					float(row[2]),
					float(row[1]),
					float(row[3])
					)
				)
		# get user home, work, study locations
		with open(config.input_survey_responses_file, newline='') as f:
                        #TODO should get passed along like the raw data
			reader = csv.DictReader(f)
			for row in reader:
				if row['uuid'] != user_id:
					continue
				if row['location_home_lat'] != '':
					self.home = Location(row['location_home_lon'],row['location_home_lat'])
				if row['location_work_lat'] != '':
					self.work = Location(row['location_work_lon'],row['location_work_lat'])
				if row['location_study_lat'] != '':
					self.school = Location(row['location_study_lon'],row['location_study_lat'])
		# sort the list by time
		self.points.sort( key=lambda x: x.epoch )
		# measure to and from neighbors
		all_indices = [ i for i,p in enumerate(self.points) ]
		self.observe_neighbors( all_indices )

	def flush(self): #TODO bottlenecked
		"""After everything is finished write all the output from this trace. 
			All writing to files should be done here if possible. Any data that 
			needs to ultimately find it's way here should be stored as a property.
			All Trace's call this at the end, and the files are initialized in main
			so we only append rows here."""
		# write potential activity locations to file
		with open(config.output_locations_file, "a") as f:
			for location in self.locations:
				f.write( "{},{},{},{},{},{}\n".format(
					self.id,				# user_id
					location.id,		# location_id
					location.longitude, 
					location.latitude, 
					location.name,		# description 
					location.visited 	# whether it was used or not
				) )
		# write episodes file
		with open(config.output_episodes_file, "a") as f:
			for i, episode in enumerate(self.episodes):
				f.write( "{},{},{},{},{},{}\n".format(
					self.id,		# user_id
					i,				# activity sequence
					episode.location_id,	# location_id
					'',						# mode
					episode.unknown,	# unknown
					episode.start		# start_time
				) )
		# write preliminary points file
		# 'user_id,lon,lat,removed,interpolated,state'
		with open(config.output_points_file,'a') as f:
			for point in self.all_interpolated_points:
				f.write( "{},{},{},{},{},{}\n".format(
					self.id,		# user_id
					point.longitude,
					point.latitude,
					'',	# removed
					'',	# interpolated
					''		# state
				) )
		# output day summary file for Steve
		days = self.get_days()
		with open(config.output_days_file,'a') as f:
			for date in days:
				#print( days[date] )
				f.write( "{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
					self.id,		# user_id
					date,
					date.weekday(),
					sum(days[date]['total']),	# total minutes
					len(days[date]['travel']),	# trip count
					sum(days[date]['travel']),	# travel time
					sum(days[date]['unknown']),# unknown time
					sum(days[date]['home']),	# home_time
					sum(days[date]['work']),	# work_time
					sum(days[date]['school']),	# school_time
					len(days[date]['home']),	# home count
					len(days[date]['work']),	# work count
					len(days[date]['school'])	# school count
				) )
				
	def get_days(self):
		"""Repartition episodes into day units."""
		day_offset = timedelta(hours=3)
		days = {}
		# for each episode except the last (always compare to next)
		for i in range(0,len(self.episodes)-1):
			# what date(s) did this occur on?
			start_date = ( self.episodes[i].start - day_offset ).date()
			end_date = ( self.episodes[i+1].start - day_offset ).date()
			date_range = (end_date-start_date).days
			# for each date touched by each activity
			for offset in range(0,date_range+1):
				date = start_date + timedelta(days=offset)
				# initialize the date the first time we see it
				if date not in days:
					days[date] = {
						'home': [], 'work': [], 'school': [], 'other': [],
						'travel': [], 'unknown': [], 'total': []
					}
				# how much of this activity occured on this date?
				# first limit start_time to start of this day
				slicer1 = datetime.combine(date, datetime.min.time()) + day_offset
				st1 = slicer1 if self.episodes[i].start < slicer1 else self.episodes[i].start
				# now limit end time to end of this day
				slicer2 = datetime.combine(date + timedelta(days=1), datetime.min.time()) + day_offset
				st2 = slicer2 if self.episodes[i+1].start > slicer2 else self.episodes[i+1].start
				# get duration in minutes from timedelta obj
				duration = (st2-st1).total_seconds() / 60
				# add activity duration to the total for this date
				days[date]['total'].append( duration )
				# now also add the duration to the appropriate category
				if self.episodes[i].e_type == 'trip':
					days[date]['travel'].append( duration )
				elif self.episodes[i].e_type == 'unknown':
					days[date]['unknown'].append( duration )
				elif self.episodes[i].e_type == 'activity':
					if self.episodes[i].a_type:
						# activity location in ['home','work','school']
						days[date][self.episodes[i].a_type].append( duration )
					else: # activity location not known
						days[date]['other'].append( duration )
		return days

	def interpolate_segment(self, segment):
		"""Takes a known subset (a list of ordered points) and interpolates
			linearly between them such that the distance between new points
			is never greater than a value specified in config, e.g. 30m."""
		new_points = []
		# For each point bu the last
		for i in range( 0, len(segment)-1 ):
			pair_int = segment[i].pair_interpolation(segment[i+1])
			new_points.extend(pair_int)
		new_points.append(segment[-1])
		return new_points
                
	def make_known_subsets(self):
		"""Partition the trace points into sets for which we're confident 
			we don't have substantial missing data. That is, exclude segments 
			where it seems like we have no data, but substantial movement; for 
			which trip and activity reconstruction would be impossible.
			TODO: Eventually we will need some much stricter checking here 
			and eventually an explicit check foro subway trips."""
		known_segments = []
		segment = [ self.points[0] ]
		# iterate over all points (except the first). Test each point to see 
		# whether we add it to the current segment or the one after.
		for i in range(1,len(self.points)):
			if ( 
				# distance over 1 km?
				distance( self.points[i], self.points[i-1] ) > 1000 and
				# time gap > 2 hours?
				self.points[i].epoch - self.points[i-1].epoch > 1*3600
			):
				# append point to next segment
				known_segments.append( segment )
				segment = [ self.points[i] ]
			else:
				# add this point to the current segment
				segment.append( self.points[i] )
		if len(segment) > 1:
			known_segments.append( segment )
		# check these segments for plausibility and append to the global property
		for segment in known_segments:
			if (
				# at least one point
				len(segment) > 1 and
				# sufficient time between last and first points
				segment[-1].epoch - segment[0].epoch > 3600
			): # mininum time length of segment?
				self.known_subsets.append(segment)
		print( '\t',len(self.known_subsets)-1,'gap(s) found in data' )
		
	def get_activity_locations(self):
		"""Get activity locations for this trace. ( Create inputs for a KDE
			function and find peaks in the surface. )"""
		for subset in self.known_subsets:
			# interpolate the subset and weight the points
			interpolated_subset = self.interpolate_segment(subset)
			self.known_subsets_interpolated.append( interpolated_subset )
			self.weight_points( interpolated_subset )
		# get all (real & interpolated) points in one big list
		self.all_interpolated_points = [ p for s in self.known_subsets_interpolated for p in s ]
		if len(self.all_interpolated_points) > 75000:
			raise Exception('Too many points for efficient KDE')
		# format as vectors for KDE function
		Xvector = [ p.x for p in self.all_interpolated_points ]
		Yvector = [ p.y for p in self.all_interpolated_points ]
		Wvector = [ p.weight for p in self.all_interpolated_points ]
		# run the KDE
		estimates, locations = kde(Xvector,Yvector,Wvector)
		# determine average GPS accuracy value for this user
		# (sqrt of the mean variance)
		mean_accuracy = sqrt(
			sum( [p.accuracy**2 for p in self.points] )
			/ len(self.points)
		)
		# estimate peak threshold value
		threshold = min_peak(
			mean_accuracy,		# mean sd of GPS accuracy for user
			Sum(Wvector),		# total seconds entering KDE
		)
		# Find peaks in the density surface
		locations = self.find_peaks(estimates,locations,threshold)
		# store the result
		self.locations.extend( locations )
		return self.locations

	def identify_locations(self):
		"""Identify locations with user-provided home, work, school locations if 
			possible."""
		# TODO this function was written in a hurry and should be made more robust
		if self.home:
			for location in self.locations:
				if distance( location, self.home ) <= 150: # meters 
					location.identify('home')
		if self.work:
			for location in self.locations:
				if distance( location, self.work ) <= 150: # meters 
					location.identify('work')
		if self.school:
			for location in self.locations:
				if distance( location, self.school ) <= 150: # meters 
					location.identify('school')

	def break_trips(self):
		"""Use a Hidden Markov Model to classify all points as deriving from 
			either time spent travelling or time spent at one of the potential 
			activity locations. Allocate time to these sequences of activities 
			accordingly."""
		for point in self.all_interpolated_points:
			# get first-pass emission probabilities from all locations 
			dists = [ distance(loc,point) for loc in self.locations ]
			dists = [ d - 100 for d in dists ]
			dists = [ 0 if d < 0 else d for d in dists ]
			point.emit_p = [ gaussian(d,100) for d in dists ]
			# standardize to one if necessary
			if sum(point.emit_p) > 1:
				point.emit_p = [ p / sum(point.emit_p) for p in point.emit_p ]
			# prepend travel probability as the difference from 1 if there is any
			point.emit_p = [1 - sum(point.emit_p)] + point.emit_p
		# make a list of starting probabilities (50/50 start travelling, start stationary)
		start_p = [0.5]+[(0.5/len(self.locations))]*len(self.locations)
		# get list of potential state indices
		# 0 is travel, others are then +1 from their list location
		states = range(0,len(self.locations)+1)
		# list of locations that actually get used
		used_locations = set()
		# simple transition probability matrix e.g.:
		#     0   1   2	
		# 0  .8  .1  .1
		# 1  .2  .8  .0
		# 2  .2  .0  .8
		trans_p = []
		for s0 in states:
			trans_p.append([])
			for s1 in states:
				if s0 + s1 == 0: # travel -> travel
					trans_p[s0].append( 0.8 )
				elif s0 == 0: # travel -> place
					trans_p[s0].append( 0.2 / len(self.locations) )
				elif s1 == 0: # place -> travel
					trans_p[s0].append( 0.2 )
				elif s0 == s1: # place -> same place
					trans_p[s0].append( 0.8 )
				else: # place -> place (no travel)
					trans_p[s0].append( 0.0 ) 
		# run the viterbi algorithm on each known subset 
		for points in self.known_subsets_interpolated:
			# VITERBI ALGORITHM
			V = [{}]
			path = {}
			for state in states:
				# Initialize base cases (t == 0)
				V[0][state] = start_p[state] * points[0].emit_p[state]
				path[state] = [state]
			for t in range(1,len(points)):	# Run Viterbi for t > 0
				V.append({})
				newpath = {}
				for s1 in states:
					(prob, state) = max(
						[ ( V[t-1][s0] * trans_p[s0][s1] * points[t].emit_p[s1], s0 ) for s0 in states ]
					)
					V[t][s1] = prob
					newpath[s1] = path[state] + [s1]
				path = newpath	# Don't need to remember the old paths
			(prob, final_state) = max( [ (V[len(points)-1][s], s) for s in states ] )
			# get the optimal sequence of states
			state_path = path[final_state]
			# note which locations have been used TODO move this elsewhere
			for visited_id in set([s-1 for s in state_path if s != 0]):
				self.locations[visited_id].visited = True
			# now record the first activity episode
			self.episodes.append( Episode( 
				points[0].time, # start_time
				# location, if any
				None if state_path[0] == 0 else self.locations[state_path[0]-1],
				False # unknown_time
			) )
			prev_state = state_path[0]
			start_time = points[0].epoch
			# for each point after the first
			for i in range(1,len(state_path)):
				# look for state changes
				if prev_state != state_path[i]:
					# assume the transition happens halfway between points
					transition_time = points[i-1].time + (points[i].time-points[i-1].time)/2
					# output start of this new state
					# now record the first activity episode
					self.episodes.append( Episode( 
						points[i].time, # start_time
						# location, if any
						None if state_path[i] == 0 else self.locations[state_path[i]-1],
						False # unknown_time
					) )
					prev_state = state_path[i]
			# now record the first activity episode
			self.episodes.append( Episode( 
				points[-1].time,	# start_time				
				None,	# no location,
				True	# unknown time ends every segment
			) )
		print( '\tFound',len(self.episodes),'activities/trips' )

	def find_peaks(self,estimates,locations,threshold):
		"""PDF was estimated at a selection of points, which are here given as a 
			list of P values (estimates) and a list of (x,y) locations. The idea 
			is to toss out any values below the threshold and identify spatial 
			clusters among those that remain. In each such cluster, the highest 
			value is the activity location."""
		assert len(estimates) == len(locations)
		# drop values below the threshold
		locations = [ (x,y) for (x,y),est in zip(locations,estimates) if est >= threshold ]
		estimates = [ est for est in estimates if est >= threshold ]
		assert len(estimates) == len(locations)
		print('\tClustering',len(estimates),'points above',threshold,'threshold')
		if len(estimates) > 5000:
			raise Exception('distance matrix will be too large')
		# now calculate a distance-based connectivity matrix between all these points
		neighbs = []
		for i,(x1,y1) in enumerate(locations):
			neighbs.append([])
			for j,(x2,y2) in enumerate(locations):
				# use euclidian distance since this is already projected
				connection = sqrt((x1-x2)**2 + (y1-y2)**2) < config.cluster_distance
				neighbs[i].append( connection )
		# clusters will be a list of disjoint sets
		clusters = []
		# now for each point, check for cluster membership and add and merge clusters
		for i in range(0,len(neighbs)):
			# get a set of neighbor indices within distance, 
			# including this point itself ( dist = 0 )
			neighbors = set( [ j for j,n in enumerate(neighbs[i]) if n ] )
			# create list to keep track of clusters this set belongs to
			member_cluster_indices = []
			for i,cluster in enumerate(clusters):
				# check each cluster for overlap with this set
				if not cluster.isdisjoint( neighbors ):
					member_cluster_indices.append(i)
			if len(member_cluster_indices) == 0:
				#  we have no overlap, so this becomes a new cluster
				clusters.append(neighbors)
			elif len(member_cluster_indices) > 0:
				# we have one or more matching clusters
				for i in reversed(member_cluster_indices):
					# add everyhting together in a new cluster and
					# drop off the old clusters which are now merged in the new one
					neighbors = neighbors | clusters.pop(i)
					# add the new cluster
				clusters.append(neighbors)
		print( '\tFound',len(clusters),'activity locations' )
		potential_activity_locations = []
		# find the maximum estimate and a location with that value
		for cluster_index, cluster in enumerate(clusters):
			peak_height = max( [estimates[i] for i in cluster] )
			for i in cluster:
				# if this is the highest point
				if estimates[i] == peak_height:
					x,y = locations[i]
					lon,lat = unproject(x,y)
					# create a location and append to the list
					location = Location(lon,lat,cluster_index)
					potential_activity_locations.append( location )
					break
		return potential_activity_locations

	def weight_points(self,segment):
		"""Assign time-based weights to a series of sequential points.
			Values are in seconds, and are split between neighboring points.
			e.g.  |--p2-time---|
			p1----|----p2------|------p3
		"""
		assert len(segment) > 1
		# set weights of middle points
		for i in range(1, len(segment)-1):
			w1 = (segment[i].time - segment[i-1].time).total_seconds() / 2
			w2 = (segment[i+1].time - segment[i].time).total_seconds() / 2
			segment[i].add_weight(w1 + w2)
		# set weights of first and last points
		segment[0].add_weight((segment[1].time - segment[0].time).total_seconds() / 2)
		segment[-1].add_weight((segment[-1].time - segment[-2].time).total_seconds() / 2)

	def time_at_loc(self, locations, inted):
		for p in inted:
			for l in locations:
				if not p.far_from(l):
					l.time_at += p.weight

	#
	# CLEANING FUNCTIONS BELOW
	# 

	def pop_point(self, key):
		"""Pop a point off the current list and add it to the discard bin.
			Then update it's former neighbors in the list."""
		# pop and append 
		point = self.points.pop(key)
		self.discarded_points.append(point)
		# now using the (former) key, update the (former) neighbors 
		i_ante = key-1
		i_post = key # has already shifted over
		self.observe_neighbors( [i_ante,i_post] )

	def find_duplicates(self):
		"""find any repeated point locations"""
		# dictionary of unique locations, with a coordinate string as key
		# containing a list of points with that exact location
		locations = {}
		for point in self.points:
			key = str(point.geom)
			if key not in locations:
				locations[key] = [ point ]
			else:
				locations[key].append( point )
		for k in locations.keys():
			if len(locations[k]) > 1:
				print( k,locations[k] )
		
	def observe_neighbors(self,indices=[]):
		"""Get angle and distance measurements to and from adjacent points.
			Store in the point object. Operates on points, the inices of which 
			are given in a list referring to the current point list."""
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
				if i_ante <= 0: break
				i_ante -= 1
			# then forward to the next different point
			i_post = i+1
			while point.geom == self.points[i_post].geom:
				if i_post >= len(self.points)-1: break
				i_post += 1
			# distance to previous point
			point.d_ante = distance(point,self.points[i_ante])
			# distance to next point
			point.d_post = distance(point,self.points[i_post])
			# if either distance is zero, this is an end point and has no angle
			if point.d_ante * point.d_post == 0:
				point.angle = 180 
			else:
				# calculate the inner angle
				point.angle = inner_angle_sphere( self.points[i_ante], point, self.points[i_post] )
			# is point is identical with both neighbors?
			if (
				point.geom == self.points[i-1].geom and 
				point.geom == self.points[i+1].geom 
			): 
				point.inter = True

	def remove_sequential_duplicates(self):
		"""remove points where both neighbors have identical locations"""
		to_remove = []
		for i, point in enumerate(self.points):
			# if there is no distance between this and neighboring points
			if point.inter:
				to_remove.append(i)
		# remove the points from the main list to the recycling bin
		for i in reversed(to_remove):
			self.pop_point(i)
		print( '\t',len(to_remove),'points removed as duplicate' )

	def remove_known_error(self,error_limit):
		"""remove points reporting positional error beyond a given limit (meters)"""
		to_remove = []
		for i, point in enumerate(self.points):
			if point.accuracy > error_limit:
				to_remove.append(i)
		# remove the points from the main list to the recycling bin
		for i in reversed(to_remove):
			self.pop_point(i)
		print( '\t',len(to_remove),'points removed as high stated error' )

	def remove_positional_error(self):
		"""use angle/distance based cleaning rules"""
		i = self.find_error_index()
		count = 0
		while i:
			self.pop_point(i)
			i = self.find_error_index()
			count += 1
		print( '\t',count, 'points removed by positional cleaning' )

	def find_error_index(self):
		"""returns the index of the craziest point"""
		# check first for angle == 0, as these are all obviously crazy
		for i, point in enumerate(self.points):
			if point.angle == 0:
				return i
		# check for very low angle, moderately high distance
		# and return the craziest first
		errors = {}
		for i, point in enumerate(self.points):
			# all but first and last
			if i == 0 or i == len(self.points)-1: 
				continue
			nearest = min( point.d_ante, point.d_post )
			if point.angle < 15 and nearest > 100:				
				error_index = 1 / ( point.angle / nearest )
				errors[error_index] = i
			if len(errors.keys()) > 0:
				return errors[max(errors.keys())]
		return False
