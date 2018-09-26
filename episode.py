class Episode(object):
	"""
	An activity episode, which may be either travel or it may be time spent at
	some location. If there is an associated location it should be the latter.
	"""

	def __init__(self, start_time, location=None, is_unknown_time=False):
		# set initially:
		self.start = start_time        # datetime object
		self.location = location       # location object ref
		self.unknown = is_unknown_time # boolean
		self.next_episode = None       # the episode following this one if any

	def __lt__(self, other):
		"""for sorting by time"""
		return self.start < other.start

	def __str__(self):
		return 'Episode starting at '+str(self.start)

	@property
	def location_id(self):
		"""Return the id of the location if any."""
		return self.location.id if self.location else None

	@property
	def type(self):
		"""
		What kind of episode is this?
		must be in in {'unknown','travel','home','work','school','other'}
		"""
		if self.unknown:
			return 'unknown'
		elif not self.location:
			return 'travel'
		elif len(self.location.name) < 1:
			return 'other'
		else: # has a named location
			if 'home' in self.location.name: return 'home'
			if 'work' in self.location.name: return 'work'
			if 'school' in self.location.name: return 'school'
			assert False # must return one of the above options

	@property
	def duration(self):
		"""Return the duration in seconds if known."""
		if self.start and self.next_episode:
			return (self.start - self.next_episode.start).total_seconds()
		else:
			return 0

	def link_subsequent_episode(self,episode):
		"""Provide the episode following this episode so we know when it ends."""
		self.next_episode = episode

