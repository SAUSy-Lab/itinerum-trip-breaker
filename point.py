from misc_funcs import parse_ts, ts_str, distance, project, unproject
import datetime, math, config

class Point(object):
	"""A location/time point ie GPS point."""

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
		return distance(self, next) > (config.cluster_distance / 2) #and (next.time - self.time).seconds > 7200 
        
	def __repr__(self):
		return str(project(self.longitude, self.latitude))

	def __eq__(self, other):
		return self.latitude == other.latitude and self.longitude == other.longitude

	def copy(self):
		return Point(self.ts, self.longitude, self.latitude, self.accuracy, self.other_fields)

	def pair_interpolation(self, other, sample):
		new_points = [self.copy()]
		dist = distance(self, other)
		if dist > sample:
			time_dif = (other.time - self.time).seconds
			n_segs = math.ceil(dist / sample)
			seg_len = dist // n_segs
			x1, y1 = project(self.longitude, self.latitude)
			x2, y2 = project(other.longitude, other.latitude)
			dx, dy = (x2 - x1)/n_segs, (y2 - y1)/n_segs
			dt = time_dif / n_segs
			for np in range(1, n_segs):
				x0, y0 = x1 + np*dx, y1 + np*dy
				lng, lat = unproject(x0, y0)
				tstamp = self.time + datetime.timedelta(seconds=dt*np)
				ts = ts_str(tstamp, self.ts[-5:])
				new_points.append(Point(ts, lng, lat, self.accuracy, None))
		return new_points

	def add_weight(self, w):
		self.weight = w

	def close_to(self, location):
		pass
