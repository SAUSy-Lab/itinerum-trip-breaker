from misc_funcs import parse_ts, ts_str, distance, project, unproject
from datetime import timedelta
from math import ceil
import config


class Point(object):
	"""
	A location/time point ie GPS point.
	"""

	def __init__(self, timestamp, longitude, latitude, accuracy_meters):
		# set initially:
		self.accuracy = accuracy_meters
		self.latitude = latitude
		self.longitude = longitude
		self.X = None  # do not access this directly
		self.Y = None  # do not access this directly
		# a string representation of a timestampt
		self.ts = timestamp
		# datetime representation of the same timestamp
		self.time, self.tz = parse_ts(timestamp)
		# these  get set later... just defining them here
		self.d_ante = None  # distance to previous point
		self.d_post = None  # distance to next point
		self.angle = None  # angle between this and adjacent points
		self.inter = False  # point shares location with both neighbors?
		self.error_index = 0  # measure used during point cleaning
		self.weight = 0  # time-based weight for KDE function
		self.emit_p = []  # emission probabilities for set of travel + locations
		self.state = None  # HMM state classification
		self.kde_p = None  # estimated PDF at this point
		# for diagnostic output
		self.discarded = False  # will be true if point tossed in cleaning
		self.synthetic = False  # was this point generated e.g. by interpolation?

	@property
	def geom(self):
		"""
		Used basically to check location uniqueness.
		"""
		return (self.latitude, self.longitude)

	@property
	def epoch(self):
		"""
		Return the time in seconds since the epoch.
		"""
		return self.time.timestamp()

	@property
	def x(self):
		"""
		Return the projected X value.
		"""
		if not self.X:
			self.project()
		return self.X

	@property
	def y(self):
		"""
		Return the projected Y value.
		"""
		if not self.Y:
			self.project()
		return self.Y

	def project(self):
		"""
		Set projected x,y values from lon,lat.
		"""
		self.X, self.Y = project(self.longitude, self.latitude)

	def __repr__(self):
		return str(project(self.longitude, self.latitude))

	def __eq__(self, other):
		return self.latitude == other.latitude and self.longitude == other.longitude

	def copy(self):
		return Point(self.ts, self.longitude, self.latitude, self.accuracy)

	def pair_interpolation(self, other_point):
		"""
		Given this and one other Point object, attempts to supply a list of
		interpolated points between the two such that gaps between the points
		are never greater than config.interpolation_distance.
		"""
		new_points = [self.copy()]  # TODO Why is this copied?
		dist = distance(self, other_point)
		if dist > config.interpolation_distance:
			time_dif = (other_point.time - self.time).seconds
			n_segs = ceil(dist / config.interpolation_distance)
			seg_len = dist // n_segs
			x1, y1 = project(self.longitude, self.latitude)
			x2, y2 = project(other_point.longitude, other_point.latitude)
			dx, dy = (x2 - x1)/n_segs, (y2 - y1)/n_segs
			dt = time_dif / n_segs
			for np in range(1, n_segs):
				x0, y0 = x1 + np*dx, y1 + np*dy
				lng, lat = unproject(x0, y0)
				tstamp = self.time + timedelta(seconds=dt*np)
				ts = ts_str(tstamp, self.ts[-5:])
				new_point = Point(ts, lng, lat, self.accuracy)
				new_point.synthetic = True
				new_points.append(new_point)
		return new_points

	def weight_decimal(self, param):
		assert param > 0
		return (-1) / (param + 2) + 1

	def add_weight(self, weight):
		"""
		Assigns time-based weight value.
		"""
		assert weight >= 0  # we may want negative weights eventually
		self.weight = weight
