class Episode(object):
	"""
	An activity episode, which may be either travel or it may be time spent at
	some location. If there is an associated location it should be the latter.
	"""

	def __init__(self, start_time, location=None, is_unknown_time=False):
		# set initially:
		self.start = start_time         # datetime object
		self.location = location        # location object ref
		self.unknown = is_unknown_time  # boolean
		self.end = None                 # datetime object
		self.duration = None            # timedelta object

	def __lt__(self, other):
		"""for sorting by time"""
		return self.start < other.start

	@property
	def location_id(self):
		"""Return the id of the location if any."""
		return self.location.id if self.location else None

	@property
	def type(self):
		"""
		What kind of episode is this?
		must be in in ['unknown','travel','home','work','school','other']
		"""
		if self.unknown:
			return 'unknown'
		elif not self.location:
			return 'travel'
		elif self.location.name in ['home', 'work', 'school']:
			return self.location.name
		else:
			return 'other'
