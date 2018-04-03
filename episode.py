class Episode(object):
	"""An activity episode. May be stationary or travel. May have associated 
		locations and definitely has a start and end time."""

	def __init__(self, start_time, location, unknown ):
		# set initially:
		self.start = start_time		# datetime object
		self.location = location	# location ref
		self.unknown = unknown		# boolean
		self.end = None				# datetime object
		self.duration = None			# timedelta object

	@property
	def location_id(self):
		"""Return the id of the location if any, otherwise None."""
		return self.location.id if self.location else None
