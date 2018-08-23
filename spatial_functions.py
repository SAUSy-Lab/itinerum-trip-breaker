# functions dealing with geometry, space, projectsions, etc

from geopy.distance import great_circle
from pyproj import Proj, transform
from math import radians, sin, cos, atan2, degrees


def project(longitude, latitude, projection_string='epsg:3347'):
	"""
	Project lat-lon values. Default of 3347 is StatsCan Lambert.
	Units in meters.
	"""
	inProj = Proj(init='epsg:4326')
	outProj = Proj(init=projection_string)
	x, y = transform(inProj, outProj, longitude, latitude)
	return x, y


def unproject(x, y, from_projection_string='epsg:3347'):
	"""
	Unproject to lat-lon values. Default of 3347 is StatsCan Lambert.
	"""
	inProj = Proj(init=from_projection_string)
	outProj = Proj(init='epsg:4326')
	longitude, latitude = transform(inProj, outProj, x, y)
	return longitude, latitude


def distance(point1, point2, euclid=False):
	"""
	Gives the great circle distance between two point objects.
	Returns meters.
	"""
	if euclid:
		return sqrt((point1.X-point2.X)**2 + (point1.Y-point2.Y)**2)
	else:
		# format the inputs
		p1 = (point1.latitude, point1.longitude)
		p2 = (point2.latitude, point2.longitude)
		return great_circle(p1, p2).meters


def inner_angle_sphere(point1, point2, point3):
	"""
	Given three point objects, calculate      p1
	the smaller of the two angles formed by    \    a
	the sequence using great circles.           \
	Returns degrees.                             p2------p3
	Latitude and Longitude attributes must be available (in degrees).
	"""
	# be sure there IS an angle to measure
	assert point1 != point2 and point2 != point3
	# first compass bearing from 2 -> 1
	lat1 = radians(point1.latitude)
	lat2 = radians(point2.latitude)
	diffLong = radians(point1.longitude - point2.longitude)
	x = sin(diffLong) * cos(lat1)
	y = cos(lat2) * sin(lat1) - (sin(lat2) *
		cos(lat1) * cos(diffLong))
	bearing1 = (degrees(atan2(x, y))+360) % 360
	# second compass bearing from 2 -> 3
	lat2 = radians(point2.latitude)
	lat3 = radians(point3.latitude)
	diffLong = radians(point3.longitude - point2.longitude)
	x = sin(diffLong) * cos(lat3)
	y = cos(lat2) * sin(lat3) -\
		(sin(lat2) * cos(lat3) * cos(diffLong))
	bearing2 = (degrees(atan2(x, y))+360) % 360
	# we want the smaller of the two angles
	degree_difference = min(abs(bearing1-bearing2),
				(360 - abs(bearing1-bearing2)))
	assert degree_difference <= 180
	return degree_difference
