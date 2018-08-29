# functions dealing with geometry, space, projectsions, etc

from geopy.distance import great_circle
from math import radians, sin, cos, atan2, degrees


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
