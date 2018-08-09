from misc_funcs import ts_str, distance, project, unproject
from datetime import timedelta, datetime
from pytz import timezone
from math import ceil
import config, re

localTime = timezone(config.local_timezone)

class Point(object):
	"""
	A location/time point ie GPS point.
	"""

	def __init__(self, time, longitude, latitude, accuracy_meters):
		# set initially:
		self.latitude = float(latitude)
		self.longitude = float(longitude)
		self.accuracy = float(accuracy_meters)
		# timezone aware datetime object
		self.time = localTime.localize( datetime.fromtimestamp(float(time)) )
		# these get set later; just defining them here for clarity
		self.X = None           # do not access this directly
		self.Y = None           # do not access this directly
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
		self.known_subset = None # known subset to which this belongs if any
		self.discarded = False   # will be true if point tossed in cleaning
		self.synthetic = False   # was this point synthesized e.g. by interpolation?

	@property
	def geom(self):
		"""
		Used basically to check location uniqueness.
		"""
		return (self.latitude, self.longitude)

	@property
	def unix_time(self):
		"""
		Return the Unix time (seconds since the UTC epoch).
		"""
		return self.time.timestamp()

	@property
	def human_timestamp(self):
		"""
		Returns a human readable timestamp string.
		"""
		return self.time.isoformat(' ')

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

	def __hash__(self):
		return id(self)

	def copy(self):
		return Point(self.unix_time, self.longitude, self.latitude, self.accuracy)


	def weight_decimal(self, param):
		"""
		Retuns a value varying from 0.5 to 1 as param varies from 0 to infinity.
		Allows weight to be distributed between points according to
		a exponential function in the ratio of their distance and 
		time difference
		"""
		return (-1) / (param + 2) + 1

	def add_weight(self, weight):
		"""
		Assigns time-based weight value.
		"""
		self.weight += weight
