import config
import osmium
from points import Point


class Map(object):
	"""
	Contains any data about the outside world, i.e. not from the survey itself. 
	This is currently meant to hold only OSM data for the purposes of identifying 
	subway trips but could later contain functions like map-matching of trips to 
	the street network etc. 
	"""
	def __init__(self):
		self.dh = DataHandler()
		self.dh.apply_file('inputs/subway_stations.osm',locations=True)
		self.subway_stations = self.dh.subway_stations

	def near_subway(self,point):
		for st_uid, st in self.subway_stations.items():
			if st.distance(point) < 200: # meters
				return True
		return False


class DataHandler(osmium.SimpleHandler):
	def __init__(self):
		osmium.SimpleHandler.__init__(self)
		self.subway_stations = {}
		self.ways = {}

	def node(self, n):
		if 'railway' in n.tags and n.tags['railway'] == 'station':
			if 'station' in n.tags and n.tags['station'] == 'subway':
				# coordinates are given as large integers
				lon = int(n.location.x) / 10**7
				lat = int(n.location.y) / 10**7
				self.subway_stations[n.id] = Point(lon,lat)

	def way(self, w):
		return
		# has w.id, w.nodes, w.tags
		#nodes = []
		#for n in w.nodes:
		#	nodes.append(n.ref)
		#self.ways[w.id] = {'nodelist':nodes}

	def relation(self, r):
		# has r.id, r.members, r.tags
		pass


# initialize the map with data when the module is called
# this prevents each call from reloading OSM data which will be shared
tube_map = Map()

