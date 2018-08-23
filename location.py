class Location(object):
	"""
	An activity location, defined as a point, possibly with a name
	"""

	def __init__(self, longitude, latitude, id_num=None):
		# set initially:
		self.latitude = float(latitude)
		self.longitude = float(longitude)
		self.id = id_num
		# total time spent at this location in seconds
		self.total_time_at = 0
		# whether or not there is an actual activity recored here
		self.visited = False
		self.name = ''

	def add_time(self, seconds):
		"""
		Add seconds spent at this location.
		"""
		assert seconds >= 0
		self.total_time_at += seconds

	def identify(self, name):
		"""
		Name this location e.g. ['home','work','school'].
		"""
		self.name = name

	def __hash__(self):
		return hash(hex(id(self)))

	def __str__(self):
		return "{}, {}".format(self.latitude, self.longitude)

	def __eq__(self, other):
		return (type(other) == type(self) and
			self.latitude == other.latitude and
			self.longitude == other.longitude)
