class ActivityLocation(object):
	"""A location point"""

	def __init__(self, longitude, latitude, id_num):
		# set initially:
		self.latitude = latitude
		self.longitude = longitude
		self.id = id_num
		self.x = None
		self.y = None
		self.time_at = 0
