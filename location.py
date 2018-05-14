class Location(object):
	"""A location point"""

	def __init__(self, longitude, latitude, id_num = None):
		# set initially:
		self.latitude = float(latitude)
		self.longitude = float(longitude)
		self.id = id_num
		# total time spent at this location in seconds
		self.total_time_at = 0
		# whether or not there is an actual activity recored here
		self.visited = False
		self.name = ''

	def add_time(self,seconds):
		"""Add seconds spent at this location."""
		assert seconds >= 0
		self.total_time_at += seconds

	def identify(self,name):
		"""Name this location ('home','work','school')."""
		self.name = name
