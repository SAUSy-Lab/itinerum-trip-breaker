class ActivityLocation(object):
	"""a location point"""

	def __init__(self, longitude, latitude):
		# set initially:
		self.latitude = latitude
		self.longitude = longitude
		self.x = None
		self.y = None
