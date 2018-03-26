#
# This file defines functions not associated with object classes
#

import math

def min_peak(GPS_error_sd,kernel_sd,total_time,threshold_time):
	"""Estimate minimum peak height given time threshold and variance parameters.
		We assume that the volume under the total KDE PDF ~= 1. 
		-	total_time is the sum of the time weights used in the KDE.
		-	threshold_time is the minimum activity time.
		Times are given in seconds"""
	from scipy.stats import multivariate_normal
	total_variance = GPS_error_sd**2 + kernel_sd**2
	peak_height = multivariate_normal.pdf(
		[0.5,0.5],	# quantiles (center)
		[0,0],		# center 
		[total_variance,total_variance]	# covariance matrix 
	)
	# this is the peak height if we have no movement and 
	# total_time == threshold_time
	assert total_time > threshold_time
	return peak_height * (float(threshold_time) / total_time)


def kde(x_vector,y_vector,weights,bandwidth):
	"""Do weighted 2d KDE in R KS package, returning python results.
		Returns two lists: P estimates and estimate locations as x,y tuples."""
	# check the inputs
	assert len(x_vector) == len(y_vector)
	assert len(weights) == len(x_vector)

	# GEOTESTING: checking input geometry
	import csv
	with open('outputs/TESTING_kde_inputs.csv', 'w+') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quotechar='"')
		writer.writerow(['x','y','w'])
		for x,y,w in zip(x_vector,y_vector,weights):
			writer.writerow([x,y,w])

	# normalize the weights to the sample size
	if sum(weights) != len(weights):
		adjust_factor = len(weights) / float(sum(weights))
		weights = [ w * adjust_factor for w in weights ]
	# get the ks package with kde function
	from rpy2.robjects.packages import importr
	ks = importr('ks')
	# get basic R functions into Python
	from rpy2.robjects import r
	cbind = r['cbind']
	diag = r['diag']
	# R data type conversion
	from rpy2.robjects import FloatVector
	# do the KDE
	print( '\tRunning KDE on',len(x_vector),'points' )
	point_matrix = cbind( FloatVector(x_vector), FloatVector(y_vector) )
	surface = ks.kde(
		# points and evaluation points are the same
		x = point_matrix,
		eval_points = point_matrix,
		# weights
		w = FloatVector( weights ),
		# bandwidth / covariance matrix
		H = diag( FloatVector( [ bandwidth**2, bandwidth**2 ] ) )
	)
	eval_points = surface.rx2('eval.points')
	estimates = surface.rx2('estimate')
	# turn these into more pythonish objects so that the rpy2 syntax doesn't 
	# have to leave this function
	eva, est = [], []
	for i in range(1,len(weights)+1):
		# insert estimate values
		est.append(estimates.rx(i)[0])
		# insert location tuples
		eva.append( ( eval_points.rx(i,True)[0], eval_points.rx(i,True)[1] ) )
	# these are now vectors (python lists) giving estimated probabilities
	# and locations as x,y tuples

	# GEOTESTING: checking sampling geometry
	import csv
	with open('outputs/TESTING_kde-eval-points.csv', 'w+') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quotechar='"')
		writer.writerow(['x','y','estimate'])
		for p,(x,y) in zip(est,eva):
			writer.writerow([x,y,p])

	return est, eva


def project(longitude,latitude,projection_string='epsg:3347'):
	"""Project lat-lon values. Default of 3347 is StatsCan Lambert.
		Units in meters."""
	from pyproj import Proj, transform
	inProj = Proj( init = 'epsg:4326' )
	outProj = Proj( init = projection_string )
	x,y = transform( inProj, outProj, longitude, latitude )
	return x,y


def unproject(x,y,from_projection_string='epsg:3347'):
	"""Unproject to lat-lon values. Default of 3347 is StatsCan Lambert."""
	from pyproj import Proj, transform
	inProj = Proj( init = from_projection_string )
	outProj = Proj( init = 'epsg:4326' )
	longitude,latitude = transform( inProj, outProj, x, y )
	return longitude,latitude


def ts_str(ts, tz):
	"""DOCUMENTATION NEEDED"""
	mo = str(ts.month) if ts.month > 9 else "0"+str(ts.month)
	d = str(ts.day) if ts.day > 9 else "0"+str(ts.day)
	h = str(ts.hour) if ts.hour > 9 else "0"+str(ts.hour)
	mi = str(ts.minute) if ts.minute > 9  else "0"+str(ts.minute)
	s = str(ts.second) if ts.second > 9 else "0"+str(ts.second)
	return "{}-{}-{}T{}:{}:{}-{}".format(ts.year, mo, d, h, mi, s, tz)


def parse_ts(timestamp): # I need to fix this
	"""DOCUMENTATION NEEDED"""
	import datetime
	# ts = 'YYYY-MM-DDThh:mm:ss-00:00'
	year = int(timestamp[:4])
	month = int(timestamp[5:7])
	day = int(timestamp[8:10])
	hour = int(timestamp[11:13])
	minutes = int(timestamp[14:16])
	second = int(timestamp[17:19])
	tz = timestamp[20:]
	return datetime.datetime(year, month, day, hour, minutes, second), tz


def distance(point1,point2):
	"""Gives the great circle distance between two point objects.
		Returns meters."""
	# import the function...
	from geopy.distance import great_circle
	# format the inputs
	p1 = ( point1.latitude, point1.longitude )
	p2 = ( point2.latitude, point2.longitude )
	return great_circle( p1, p2 ).meters


def inner_angle_sphere(point1,point2,point3):
	"""Given three point objects, calculate      p1
		the smaller of the two angles formed by    \    a
		the sequence using great circles.           \  
		Returns degrees.                             p2------p3
		Latitude and Longitude attributes must be available (in degrees).
	"""
	# be sure there IS an angle to measure
	assert point1 != point2 and point2 != point3
	# first compass bearing from 2 -> 1
	lat1 = math.radians(point1.latitude)
	lat2 = math.radians(point2.latitude)
	diffLong = math.radians(point1.longitude - point2.longitude)
	x = math.sin(diffLong) * math.cos(lat1)
	y = math.cos(lat2) * math.sin(lat1) - (math.sin(lat2) * math.cos(lat1) * math.cos(diffLong))
	bearing1 = (math.degrees(math.atan2(x, y))+360) % 360 
	# second compass bearing from 2 -> 3
	lat2 = math.radians(point2.latitude)
	lat3 = math.radians(point3.latitude)
	diffLong = math.radians(point3.longitude - point2.longitude)
	x = math.sin(diffLong) * math.cos(lat3)
	y = math.cos(lat2) * math.sin(lat3) - (math.sin(lat2) * math.cos(lat3) * math.cos(diffLong))
	bearing2 = (math.degrees(math.atan2(x, y))+360) % 360
	# we want the smaller of the two angles
	degree_difference = min( abs(bearing1-bearing2), (360 - abs(bearing1-bearing2)) )
	assert degree_difference <= 180
	return degree_difference
