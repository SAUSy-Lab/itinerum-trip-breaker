#Meaningless change
import config, csv, math
from point import Point
from misc_funcs import distance, inner_angle_sphere, project, kde, min_peak, gaussian, ts_str

class Trace(object):
	"""A "trace", a GPS trace, is all the data associated with one itinerum user.
		It's mainly treated here as a temporal/spatial sequence of points."""

	def __init__(self,user_id):
		"""Construct the user object by pulling all data pertaining to this user.
			Identified by ID"""
		self.id = user_id		# 
		self.points = []		# time-ordered list of points
		self.discarded_points = [] # list of points removed
		# ordered list of ordered lists of points separated by unknown times
		# this will be the basic object of any trip-breaking analysis
		self.known_subsets = []
		self.known_subsets_interpolated = []
		# list of potential activity locations
		self.locations = []
		# records the number of activity records written so far
		self.activity_count = 0
		# read in all time and location data for points
		# right now only using a few fields
		with open(config.input_coordinates_file, newline='') as f:
			reader = csv.DictReader(f)
			for row in reader:
				if row['uuid'] != user_id:
					continue
				# add point to the list
				self.points.append(
					Point(
						row['timestamp'],
						float(row['longitude']),
						float(row['latitude']),
						float(row['h_accuracy'])
					)
				)
		# measure to and from neighbors
		all_indices = [ i for i,p in enumerate(self.points) ]
		self.observe_neighbors( all_indices )

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
				self.points[i].epoch - self.points[i-1].epoch > 2*3600
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
		kde_points = [ p for s in self.known_subsets_interpolated for p in s ]
		# format as vectors for KDE function
		Xvector = [ p.x for p in kde_points ]
		Yvector = [ p.y for p in kde_points ]
		Wvector = [ p.weight for p in kde_points ]
		# run the KDE
		estimates, locations = kde(Xvector,Yvector,Wvector)
		# determine average GPS accuracy value for this user
		# (sqrt of the mean variance)
		mean_accuracy = math.sqrt(
			sum( [p.accuracy**2 for p in self.points] ) 
			/ len(self.points)
		)
		# estimate peak threshold value
		threshold = min_peak(
			mean_accuracy,		# mean sd of GPS accuracy for user
			sum(Wvector),		# total seconds entering KDE
		)
		# Find peaks in the density surface
		locations = self.find_peaks(estimates,locations,threshold)
		# store the result
		self.locations.extend( locations )
		return self.locations

	def break_trips(self):
		"""Use a Hidden Markov Model to classify all points as deriving from 
			either time spent travelling or time spent at one of the potential 
			activity locations. Allocate time to these sequences of activities 
			accordingly."""
		for point in self.points:
			# get first-pass emission probabilities from all locations 
			dists = [ distance(loc,point) for loc in self.locations ]
			dists = [ d - config.cluster_distance/2 for d in dists ]
			dists = [ 0 if d < 0 else d for d in dists ]
			point.emit_p = [ gaussian(d,25) for d in dists ]
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
		for points in self.known_subsets:
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
			#print(state_path)
			# note which locations have been used
			used_locations = used_locations | set([s-1 for s in state_path if s != 0])
			# NOW WE OUTPUT THE ACTIVITIES TO FILE
			# output the first activity start
			self.write_activity( 
				'' if state_path[0] == 0 else state_path[0] - 1, # location_id 
				'', # mode 
				False, # unknown time
				ts_str(points[0].time,points[0].tz) # timestamp
			)
			prev_state = state_path[0]
			start_time = points[0].epoch
			# for each point after the first
			for i in range(1,len(state_path)):
				# look for state changes
				if prev_state != state_path[i]:
					# assume the transition happens halfway between points
					transition_time = points[i-1].time + (points[i].time-points[i-1].time)/2
					# output start of this new state
					self.write_activity( 
						'' if state_path[i] == 0 else state_path[i] - 1, # location_id 
						'', # mode 
						False, # unknown time
						ts_str(transition_time,points[i].tz) # timestamp
					)
					prev_state = state_path[i]
			# output the last (unknown) activity to cap off this segment
			self.write_activity( 
				'', # location_id 
				'', # mode 
				True, # unknown time
				ts_str(points[-1].time,points[-1].tz) # timestamp
			)
		print( '\tFound',self.activity_count,'activities/trips' )
		# write the locations, whether used or not
		for location in self.locations:
			self.write_location(location,location.id in used_locations)
		
	def write_location(self, location, used):
		"""Write a single location record to the output file."""
		f = open(config.output_locations_file, "a")
		line = "{},{},{},{},{},{}\n".format(
			self.id, 
			location.id, 
			location.longitude, 
			location.latitude, 
			'', # description 
			used
		)
		f.write(line)
		f.close()

	def write_activity( self, location_id, mode, unknown_val, start_time ):
		"""Write a single activity record to the output CSV."""
		f = open(config.output_activities_file, "a") # append to the file
		self.activity_count += 1
		line = "{},{},{},{},{},{}\n".format(
			self.id, 
			self.activity_count, 
			location_id,
			mode,
			unknown_val,
			start_time
		)
		f.write(line)
		f.close()

	def find_peaks(self,estimates,locations,threshold):
		"""PDF was estimated at a selection of points, which are here given as a 
			list of P values (estimates) and a list of (x,y) locations. The idea 
			is to toss out any values below the threshold and identify spatial 
			clusters among those that remain. In each such cluster, the highest 
			value is the activity location."""
		assert len(estimates) == len(locations)
		from math import sqrt
		from location import ActivityLocation
		from misc_funcs import unproject
		# drop values below the threshold
		locations = [ (x,y) for (x,y),est in zip(locations,estimates) if est >= threshold ]
		estimates = [ est for est in estimates if est >= threshold ]
		assert len(estimates) == len(locations)
		print('\tClustering',len(estimates),'points above',threshold,'threshold')
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
					location = ActivityLocation(lon,lat,cluster_index)
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
