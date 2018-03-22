# This is an early algorithm for cleaning itinerum output data.
# I'm mainly trying to remove obvious positional errors and drop out 
# redundant points. After this comes trip/activity partitioning. 
# 
# The key thing here is to remove points that don't appear to be based 
# on reasonably accurate *GPS*. Id est, some points come from black-box 
# systems inside the phone, and some are GPS based but not accurate. 
# There a few things we can look for to pick these out:

# 1) GPS noise has a roughly gaussian distribution, and should very 
#		rarely hit the same exact point twice, given enough precision. 
#		Yet phones often get "stuck" on a point and repeat it. 
#		This may indicate a cell-tower or some other problem. Any 
#		point that is used two or more times should come under strict 
#		scrutiny. This is not quite yet implemented here. 

# 2) Any major jump away and back again, especially if the h_error
#		Doesn't justify the distance may indicate a non-GPS signal 
#		or a bad error estimate. Away-and-back-again points are identified 
#		by the minimum distance from neighboring points and the angle formed 
#		between the three. 

# The algorithm iterates over users and advances iteratively over each 
# user's data.  

import datetime, csv, math, scratch
import rpy2

# minimum time in seconds at a location for location detection
MIN_SECS_AT_LOC = 10*60
# kernel bandwidth in meters
BANDWIDTH = 25

def inner_angle_sphere(point1,point2,point3):
	"""Given three point objects, calculate      p1
		the smaller of the two angles formed by    \    a
		the sequence using great circles.           \  
		Returns degrees.                             p2------p3
		Latitude and Longitude attributes must be available (in degrees).
	"""
	# be sure there IS an angle to measure
	assert point1 != point2 and point2 != point3
	# first compass bearing from 2 -> 1
	lat1 = math.radians(point1.latitude)
	lat2 = math.radians(point2.latitude)
	diffLong = math.radians(point1.longitude - point2.longitude)
	x = math.sin(diffLong) * math.cos(lat1)
	y = math.cos(lat2) * math.sin(lat1) - (math.sin(lat2) * math.cos(lat1) * math.cos(diffLong))
	bearing1 = (math.degrees(math.atan2(x, y))+360) % 360 
	# second compass bearing from 2 -> 3
	lat2 = math.radians(point2.latitude)
	lat3 = math.radians(point3.latitude)
	diffLong = math.radians(point3.longitude - point2.longitude)
	x = math.sin(diffLong) * math.cos(lat3)
	y = math.cos(lat2) * math.sin(lat3) - (math.sin(lat2) * math.cos(lat3) * math.cos(diffLong))
	bearing2 = (math.degrees(math.atan2(x, y))+360) % 360
	# we want the smaller of the two angles
	degree_difference = min( abs(bearing1-bearing2), (360 - abs(bearing1-bearing2)) )
	assert degree_difference <= 180
	return degree_difference


def distance(point1,point2):
	"""Gives the great circle distance between two point objects.
		Returns meters."""
	# import the function...
	from geopy.distance import great_circle
	# format the inputs
	p1 = ( point1.latitude, point1.longitude )
	p2 = ( point2.latitude, point2.longitude )
	return great_circle( p1, p2 ).meters


class point_obj(object):
	"""a location/time point"""

	def __init__(self, timestamp, longitude, latitude, accuracy_meters, other_fields):
		# set initially:
		self.accuracy = accuracy_meters
		self.latitude = latitude
		self.longitude = longitude
		self.ts = timestamp
		self.time, self.tz = parse_ts(timestamp)
		# dictionary storing other fields that will just pass through
		self.other_fields = other_fields
		# these  get set later... just defining them here
		self.d_ante = None		# distance to previous point
		self.d_post = None		# distance to next point
		self.angle = None			# angle between this and adjacent points
		self.inter = False		# point shares location with both neighbors?
		self.error_index = 0 
		self.weight = 0
		
	@property
	def geom(self):
		"""used basically to check location uniqueness"""
		return (self.latitude,self.longitude)

	def far_from(self, next):
		# next must be a Point
		return distance(self, next) > 1000 #and (next.time - self.time).seconds > 7200 
        
	def __repr__(self):
		return str(scratch.project(self.longitude, self.latitude))

	def copy(self):
		return point_obj(self.ts, self.longitude, self.latitude, self.accuracy, self.other_fields)

	def pair_interpolation(self, other, sample):
		new_points = [self.copy()]
		dist = distance(self, other)
		if dist > sample:
			time_dif = (other.time - self.time).seconds
			n_segs = math.ceil(dist / sample)
			seg_len = dist // n_segs
			x1, y1 = scratch.project(self.longitude, self.latitude)
			x2, y2 = scratch.project(other.longitude, other.latitude)
			dx, dy = (x2 - x1)/n_segs, (y2 - y1)/n_segs
			dt = time_dif / n_segs
			for np in range(1, n_segs):
				x0, y0 = x1 + np*dx, y1 + np*dy
				lng, lat = scratch.unproject(x0, y0)
				tstamp = self.time + datetime.timedelta(seconds=dt*np)
				ts = ts_str(tstamp, self.ts[-5:])
				new_points.append(point_obj(ts, lng, lat, self.accuracy, None))
		return new_points

	def add_weight(self, w):
		self.weight = w

class trace(object):
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
		with open(input_coordinates_file, newline='') as f:
			reader = csv.DictReader(f)
			for row in reader:
				if row['uuid'] != user_id:
					continue
				# add point to the list
				self.points.append(
					point_obj(
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

	def interpolate_segment(self, segment, sample=30):
		new_points = []
		for i in range(len(segment)-1):
			pair_int = segment[i].pair_interpolation(segment[i+1], sample)
			new_points.extend(pair_int)
		new_points.append(segment[-1])
		return new_points

	def make_subsets(self):
		ss = []
		cur = [self.points[0]]
		for i in range(1, len(self.points)):
			cur.append(self.points[i])
			if self.points[i-1].far_from(self.points[i]):
				ss.append(cur[:])i
				cur = []

		for known_segment in ss:
			if len(known_segment) > 1: # mininum time length of segment?
				self.subsets.append(known_segment)


	def break_trips(self):
		ml = []
		for sl in self.subsets:
			interpolated = self.interpolate_segment(sl, 30)
			weight_points(interpolated])

		# format as vectors for KDE function
		xs = [ scratch.project(p.longitude, p.latitude)[0] for p in ml]
		ys = [ scratch.project(p.longitude, p.latitude)[1] for p in ml]
		ws = [p.weight for p in ml]
		# run the KDE
		estimates, locations = scratch.kde(xs,ys,ws,BANDWIDTH)
		# get the average GPS accuracy value
		accuracy = sum([p.accuracy for p in self.points])/len(self.points)
		print('acc:',accuracy)
		# estimate peak threshold value
		threshold = scratch.min_peak(accuracy,BANDWIDTH,sum(ws),MIN_SECS_AT_LOC)
		# Find peaks in the density surface
		# currently only testing this function
		scratch.find_peaks(estimates,locations,threshold)

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

def ts_str(ts, tz):
	mo = str(ts.month) if ts.month > 9 else "0"+str(ts.month)
	d = str(ts.day) if ts.day > 9 else "0"+str(ts.day)
	h = str(ts.hour) if ts.hour > 9 else "0"+str(ts.hour)
	mi = str(ts.minute) if ts.minute > 9  else "0"+str(ts.minute)
	s = str(ts.second) if ts.second > 9 else "0"+str(ts.second)
	return "{}-{}-{}T{}:{}:{}-{}".format(ts.year, mo, d, h, mi, s, tz)


def parse_ts(timestamp): # I need to fix this
	# ts = 'YYYY-MM-DDThh:mm:ss-00:00'
	year = int(timestamp[:4])
	month = int(timestamp[5:7])
	day = int(timestamp[8:10])
	hour = int(timestamp[11:13])
	minutes = int(timestamp[14:16])
	second = int(timestamp[17:19])
	tz = timestamp[20:]
	return datetime.datetime(year, month, day, hour, minutes, second), tz

def weight_points(segment):
	for i in range(1, len(segment)-1):
		w1 = (segment[i].time - segment[i-1].time).seconds / 2
		w2 = (segment[i+1].time - segment[i].time).seconds / 2
		segment[i].add_weight(w1 + w2)
	segment[0].add_weight((segment[1].time - segment[0].time).seconds / 2)
	segment[-1].add_weight((segment[-1].time - segment[-2].time).seconds / 2)

# Standard format so we can import this module elsewhere.
if __name__ == "__main__":

	input_coordinates_file = 'inputs/coordinates.csv'
	output_coordinates_file = 'outputs/cleaned_coordinates.csv'

	user_ids = []

	# get a list of unique user_ids to iterate over
	with open(input_coordinates_file, newline='') as f:
		reader = csv.DictReader(f)
		for row in reader:
			user_ids.append( row['uuid'] )

	user_ids = list(set(user_ids)) # why? user_ids is already a list
	print( len(user_ids),'user(s) to clean' )


        # open a file for results
	with open(output_coordinates_file, 'w', newline='') as csvfile:
		fieldnames = ['user_id','latitude','longitude','h_accuracy','timestamp','use']
		fieldnames += ['v_accuracy','id','speed','point_type']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()
		# loop over users calling all the functions for each
		for user_id in user_ids:
			user = trace(user_id)
			print( len(user.points),'points at start for',user_id )
			user.remove_known_error(100)
			user.remove_sequential_duplicates()
			user.remove_positional_error()
			# this is actually necessary again after positional cleaning
			# ( some angles == 0 )
			user.remove_sequential_duplicates()

			# GEOTESTING: checking post-cleaning geometry
			import csv
			with open('outputs/TESTING_post-cleaning-points.csv', 'w+') as csvfile:
				writer = csv.writer(csvfile, delimiter=',', quotechar='"')
				writer.writerow(['latitude','longitude','timestamp'])
				for point in user.points:
					writer.writerow([point.latitude,point.longitude,point.time])

			user.make_subsets()
			user.break_trips() 

#		# now store all the points for this user
#		for point in user.points:			
#			writer.writerow(
#				{
#					**{
#						'user_id': user_id, 
#						'latitude': point.latitude,
#						'longitude': point.longitude,
#						'h_accuracy': point.accuracy,
#						'timestamp': point.time,
#						'use': True
#					},
#					**point.other_fields
#				}
#			)
#		for point in user.discarded_points:			
#			writer.writerow(
#				{
#					**{
#						'user_id': user_id, 
#						'latitude': point.latitude,
#						'longitude': point.longitude,
#						'h_accuracy': point.accuracy,
#						'timestamp': point.time,
#						'use': False
#					},
#					**point.other_fields
#				}
#			)

	

