# these are classes representing spatial points

from pyproj import Proj, transform
from geopy.distance import great_circle
from geopy.distance import great_circle
from math import radians, sin, cos, atan2, degrees
from datetime import timedelta, datetime
from pytz import timezone
import config
localTime = timezone(config.local_timezone)


class Point:
	"""Base class for spatial point types. Can be instantiated with either
	unprojected lat/lon (preferred) or x,y in the default projection"""

	def __init__( self, longitude=None, latitude=None, X=None, Y=None ):
		# set the location from either lat/lon
		if longitude and latitude:
			self.latitude = float(latitude)
			self.longitude = float(longitude)
			self.project() # sets self.X,self.Y
		# or X,Y
		elif X and Y:
			self.X = X
			self.Y = Y
			self.unproject()
		else:
			assert False # we should never be here. Not enough coordinates supplied

	@property
	def geom(self):
		"""Used basically to check location uniqueness."""
		return (self.latitude, self.longitude)

	@property
	def x(self):
		"""Return the projected X value."""
		#if not self.X:
		#	self.project()
		return self.X

	@property
	def y(self):
		"""Return the projected Y value."""
		#if not self.Y:
		#	self.project()
		return self.Y

	def distance(self, point2, euclid=False):
		"""
		Gives the great circle distance between two point objects.
		Returns meters.
		"""
		if euclid:
			return sqrt((self.x-point2.x)**2 + (self.y-point2.y)**2)
		else:
			return great_circle(self.geom, point2.geom).meters

	def project(self,projection_string='epsg:3347'):
		"""Set projected x,y values from lon,lat. 
		Default of 3347 is StatsCan Lambert, units in meters."""
		inProj = Proj(init='epsg:4326')
		outProj = Proj(init=projection_string)
		self.X, self.Y = transform(inProj, outProj, self.longitude, self.latitude)

	def unproject(self,from_projection_string='epsg:3347'):
		"""Unproject to lat-lon values. Default of 3347 is StatsCan Lambert."""
		inProj = Proj(init=from_projection_string)
		outProj = Proj(init='epsg:4326')
		self.longitude, self.latitude = transform(inProj, outProj, self.X, self.Y)

	def inner_angle_sphere(self, point1, point2):
		"""Calculate the smaller of the two angles  p1    A
		formed between three unprojected points.     \ 
		Returns degrees.                            Self----p2
		"""
		# be sure there is an angle to measure
		assert self != point1 and self != point2
		# first compass bearing from 2 -> 1
		lat1 = radians(point1.latitude)
		lat2 = radians(self.latitude)
		diffLong = radians(point1.longitude - self.longitude)
		x = sin(diffLong) * cos(lat1)
		y = cos(lat2) * sin(lat1) - (sin(lat2) *
			cos(lat1) * cos(diffLong))
		bearing1 = (degrees(atan2(x, y))+360) % 360
		# second compass bearing from 2 -> 3
		#lat2 = radians(self.latitude)
		lat3 = radians(point2.latitude)
		diffLong = radians(point2.longitude - self.longitude)
		x = sin(diffLong) * cos(lat3)
		y = cos(lat2) * sin(lat3) - (sin(lat2) * cos(lat3) * cos(diffLong))
		bearing2 = (degrees(atan2(x, y))+360) % 360
		# we want the smaller of the two angles
		degree_difference = min(abs(bearing1-bearing2),
					(360 - abs(bearing1-bearing2)))
		assert degree_difference <= 180
		return degree_difference

	def __repr__(self):
		return "Point at lat:{},lon:{}".format(self.latitude, self.longitude)

	def __hash__(self):
		return id(self)

	def copy(self):
		return Point(self.unix_time, self.longitude, self.latitude, self.accuracy)

	def __eq__(self, other):
		return (type(other) == type(self) and
			self.latitude == other.latitude and
			self.longitude == other.longitude)


class GPSpoint(Point):
	"""A space/time point ie GPS point."""

	def __init__(self,time,accuracy_meters,lon=None,lat=None,X=None,Y=None):
		Point.__init__(self,lon,lat,X,Y)
		# immediate assignments
		self.accuracy = float(accuracy_meters)
		self.time = localTime.localize(datetime.fromtimestamp(float(time)))
		# these get set later; just defining them here for clarity
		self.d_ante = None      # distance to previous point
		self.d_post = None      # distance to next point
		self.angle = None       # angle between this and adjacent points
		self.inter = False      # point shares location with both neighbors?
		self.error_index = 0    # measure used during point cleaning
		self.weight = 0         # time-based weight for KDE function
		self.emit_p = []        # emission probabilities for set of travel+locations
		self.state = None       # HMM state classification
		self.location = None    # reference to location object point is at per state
		self.kde_p = None       # estimated PDF at this point
		# for diagnostic output
		self.known_subset = None  # known subset to which this belongs if any
		self.discarded = False   # will be true if point tossed in cleaning
		self.synthetic = False   # was this point synthesized e.g. by interpolation?

	@property
	def unix_time(self):
		"""Return the Unix time (seconds since the UTC epoch)."""
		return self.time.timestamp()

	@property
	def human_timestamp(self):
		"""Returns a human readable timestamp string."""
		return self.time.isoformat(' ')

	def add_weight(self, weight):
		"""Assigns time-based weight value."""
		self.weight += weight

	def delta_t(self, other_point):
		"""Gives the absolute difference in seconds between two points."""
		return abs(self.unix_time - other_point.unix_time)

	def mps(self, other_point):
		"""Gives the speed between two points in meters per second."""
		return self.distance(other_point) / self.delta_t(other_point)


class Location(Point):
	"""An activity location, defined as a point, possibly with a name."""

	def __init__( self, longitude, latitude, id_num=None ):
		Point.__init__(self,longitude, latitude)
		# set initially:
		self.id = id_num
		# total time spent at this location in seconds
		self.total_time_at = 0
		# whether or not there is an actual activity recored here
		self.visited = False
		self.name = ''

	def add_time(self, seconds):
		"""Add seconds spent at this location."""
		assert seconds >= 0
		self.total_time_at += seconds

	def identify(self, name):
		"""Name this location e.g. ['home','work','school']."""
		self.name = name

