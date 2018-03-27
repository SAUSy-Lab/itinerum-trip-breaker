import config, csv, math
from point import Point
from misc_funcs import distance, inner_angle_sphere, project, kde, min_peak

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
		self.subsets = []
		
		other_keys = ['id','speed','v_accuracy','point_type'] # The keys you want

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
						float(row['h_accuracy']),
						dict((k, row[k]) for k in other_keys if k in row)
					)
				)
		# measure to and from neighbors
		all_indices = [ i for i,p in enumerate(self.points) ]
		self.observe_neighbors( all_indices )

	def compute_sequence(self, locations):
		"""DOCUMENTATION NEEDED"""
		sequence = []
		cur = []
		p_loc = None #previous location
		loc = None #current location
		for p in self.points:
			cur.append(p)
			found = False
			for l in locations:
				if distance(p, l) < config.cluster_distance / 2: #unique location
					p_loc = loc
					loc = l
					found = True
			if not found:
				p_loc = loc
				loc = None

			if not loc == p_loc:
				sequence.append(cur)
				cur = []
		return sequence

	def clean_sequence(self,sequence):
		"""DOCUMENTATION NEEDED"""
		pass

	def write_a_csv(self, sequence, point_to_lid, l_to_uid, filename):
		"""DOCUMENTATION NEEDED"""
		fd = open(filename, "a") #append to the file
		s_no = 1
		for event in sequence:
			mode = ""
			unknown = ""
			time = event[0].ts
			location_id = ""
			if time in point_to_lid:
				location_id = l_to_uid[(point_to_lid[time].longitude, point_to_lid[time].latitude)]
			line = "{},{},{},{},{},{}\n".format(
			        self.id, str(s_no), location_id, mode, unknown, time)
			fd.write(line)
			s_no += 1
		fd.close()

	def write_l_csv(self, locations, filename):
		"""DOCUMENTATION NEEDED"""
		fd = open(filename, "a")
		d = {}
		uid = 1
		for location in locations:
			description = ""
			line = "{},{},{},{},{}\n".format(self.id, str(uid), str(location.longitude), str(location.latitude), description)
			fd.write(line)
			if (location.longitude, location.latitude) not in d:
				d[(location.longitude, location.latitude)] = uid
			uid += 1
		fd.close()
		return d

	def interpolate_segment(self, segment, sample=30):
		"""DOCUMENTATION NEEDED"""
		new_points = []
		for i in range(len(segment)-1):
			pair_int = segment[i].pair_interpolation(segment[i+1], sample)
			new_points.extend(pair_int)
		new_points.append(segment[-1])
		return new_points

	def make_subsets(self):
		"""DOCUMENTATION NEEDED"""
		ss = []
		cur = [self.points[0]]
		for i in range(1, len(self.points)):
			cur.append(self.points[i])
			if self.points[i-1].far_from(self.points[i]):
				ss.append(cur[:])
				cur = []
		for known_segment in ss:
			if len(known_segment) > 1: # mininum time length of segment?
				self.subsets.append(known_segment)

	def break_trips(self):
		"""DOCUMENTATION NEEDED"""
		ml = []
		for sl in self.subsets:
			interpolated = self.interpolate_segment(sl, 30)
			self.weight_points(interpolated)
			ml.extend(interpolated)
		# format as vectors for KDE function
		xs = [ project(p.longitude, p.latitude)[0] for p in ml]
		ys = [ project(p.longitude, p.latitude)[1] for p in ml]
		ws = [p.weight for p in ml]
		# run the KDE
		estimates, locations = kde(xs,ys,ws,config.kernel_bandwidth)
		# determine average GPS accuracy value for this user
		# (sqrt of the mean variance)
		mean_accuracy = math.sqrt(
			sum( [p.accuracy**2 for p in self.points] ) 
			/ 
			len(self.points)
		)
		# estimate peak threshold value
		threshold = min_peak(
			mean_accuracy,		# mean sd of GPS accuracy for user
			sum(ws),				# total seconds entering KDE
		)
		# Find peaks in the density surface
		# currently only testing this function
		locations = self.find_peaks(estimates,locations,threshold)
		sequence = self.compute_sequence(locations)
		self.clean_sequence(sequence)
		ptl = self.make_ptl(locations)
		l_to_uid = self.write_l_csv(locations, config.output_locations_file)
		self.write_a_csv(sequence, ptl, l_to_uid, config.output_activities_file)

	def make_ptl(self, locations):
		"""DOCUMENTATION NEEDED"""
		d = {}
		for p in self.points:
			for l in locations:
				if distance(p, l) < config.cluster_distance / 2:
					d[p.ts] = l
		return d

	def find_peaks(self,estimates,locations,threshold):
		"""PDF was estimated at a selection of points, which are here given as a list
			of P values (estimates) and a list of (x,y) locations. The idea is to toss 
			out any values below the threshold and identify spatial clusters among 
			those that remain. In each such cluster, the highest value is the activity 
			location."""
		assert len(estimates) == len(locations)
		from math import sqrt
		from location import ActivityLocation
		from misc_funcs import unproject
		# drop values below the threshold
		locations = [ (x,y) for (x,y),est in zip(locations,estimates) if est >= threshold ]
		estimates = [ est for est in estimates if est >= threshold ]
		assert len(estimates) == len(locations)
		print('\tclustering',len(estimates),'points above',threshold,'threshold')
		# now calculate a distance-based connectivity matrix between all these points
		neighbs = []
		for i,(x1,y1) in enumerate(locations):
			neighbs.append([])
			for j,(x2,y2) in enumerate(locations):
				# use euclidian distance since this is already projected
				connection = sqrt((x1-x2)**2 + (y1-y2)**2) < config.cluster_distance
				neighbs[i].append( connection )
		print( '\thave connection matrix with',sum([len(l) for l in neighbs]),'entries' )
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

		print( '\tfound',len(clusters),'clusters with',sum([len(c) for c in clusters]),'total points' )
		potential_activity_locations = []
		# find the maximum estimate and a location with that value
		for cluster in clusters:
			peak_height = max( [estimates[i] for i in cluster] )
			for i in cluster:
				# if this is the highest point
				if estimates[i] == peak_height:
					x,y = locations[i]
					lon,lat = unproject(x,y)
					# create a location and append to the list
					location = ActivityLocation(lon,lat)
					potential_activity_locations.append( location )
					break

		# GEOTESTING: checking post-cleaning geometry
		import csv
		with open('outputs/TESTING_potential-activity-locations.csv', 'w+') as csvfile:
			writer = csv.writer(csvfile, delimiter=',', quotechar='"')
			writer.writerow(['longitude','latitude'])
			for location in potential_activity_locations:
				writer.writerow([location.longitude,location.latitude])

		return potential_activity_locations

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

	def weight_points(self,segment):
		"""DOCUMENTATION NEEDED"""
		for i in range(1, len(segment)-1):
			w1 = (segment[i].time - segment[i-1].time).seconds / 2
			w2 = (segment[i+1].time - segment[i].time).seconds / 2
			segment[i].add_weight(w1 + w2)
		segment[0].add_weight((segment[1].time - segment[0].time).seconds / 2)
		segment[-1].add_weight((segment[-1].time - segment[-2].time).seconds / 2)
