#
# This file defines functions not associated with object classes
#
import math
import config
import datetime
from geopy.distance import great_circle
from pyproj import Proj, transform
from scipy.stats import multivariate_normal
from scipy.stats import multivariate_normal


def min_peak(GPS_error_sd, total_time):
	"""
	Estimate minimum peak height given time threshold and variance parameters.
	We assume that the volume under the total KDE PDF ~= 1.
	-	total_time is the sum of the time weights used in the KDE.
	-	threshold_time is the minimum activity time.
	Times are given in seconds
	"""
	total_variance = GPS_error_sd**2 + config.kernel_bandwidth**2
	peak_height = multivariate_normal.pdf([0.5, 0.5], [0, 0],
					[total_variance, total_variance])
	# this is the peak height if we have no movement and
	# total_time == threshold_time
	assert total_time > config.minimum_activity_time
	rv = peak_height * (float(config.minimum_activity_time) / total_time)
	return rv


def kde(x_vector, y_vector, weights):
	"""
	Do weighted 2d KDE in R KS package, returning python results.
	Returns two lists: P estimates and estimate locations as x,y tuples.

	Another possible way of doing KDE (in Python) is with
	http://pysal.readthedocs.io/en/latest/users/tutorials/
	smoothing.html#non-parametric-smoothing ?
	or with http://scikit-learn.org/stable/modules/density.html ?
	"""
	from rpy2.robjects.packages import importr  # for importing packages
	from rpy2.robjects import r, FloatVector  # variable names
	# ensure all inputs are vectors of the same length
	assert len(x_vector) == len(y_vector)
	assert len(weights) == len(x_vector)
	# normalize the weights to the sample size
	if sum(weights) != len(weights):
		adjust_factor = len(weights) / float(sum(weights))
		weights = [w * adjust_factor for w in weights]
	# load the ks package (making kde function available)
	ks = importr('ks')
	# get basic R functions into Python
	cbind, diag = r['cbind'], r['diag']
	# do the KDE
	print('\tRunning KDE on', len(x_vector), 'points')
	point_matrix = cbind(FloatVector(x_vector), FloatVector(y_vector))
	bandwidth = config.kernel_bandwidth
	surface = ks.kde(
		x=point_matrix,
		eval_points=point_matrix,
		w=FloatVector(weights), 
		H=diag(FloatVector([bandwidth**2, bandwidth**2])),
		binned=False
	)
	estimates = surface.rx2('estimate')
	# turn these into more pythonish objects so that the rpy2 syntax doesn't
	# have to leave this function
	est = []
	for i in range(1, len(weights)+1):
		# insert estimate values
		est.append(estimates.rx(i)[0])
	# this is now a vector (python list)
	return est


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


def ts_str(ts, tz):
	"""
	Return a properly formatted timestamp string from a datetime
	object and a timezone string.

	inverts parse_ts
	"""
	mo = str(ts.month) if ts.month > 9 else "0"+str(ts.month)
	d = str(ts.day) if ts.day > 9 else "0"+str(ts.day)
	h = str(ts.hour) if ts.hour > 9 else "0"+str(ts.hour)
	mi = str(ts.minute) if ts.minute > 9 else "0"+str(ts.minute)
	s = str(ts.second) if ts.second > 9 else "0"+str(ts.second)
	return "{}-{}-{}T{}:{}:{}-{}".format(ts.year, mo, d, h, mi, s, tz)


def parse_ts(timestamp):
	"""
	Return a datetime object and a timezone string given a
	properly formatted timestamp string.
	Formatted according to itinerum output.

	inverts ts_str
	"""
	if timestamp == "":
		return "", ""
	# ts = 'YYYY-MM-DDThh:mm:ss-00:00'
	year = int(timestamp[:4])
	month = int(timestamp[5:7])
	day = int(timestamp[8:10])
	hour = int(timestamp[11:13])
	minutes = int(timestamp[14:16])
	second = int(timestamp[17:19])
	tz = timestamp[20:]
	return datetime.datetime(year, month, day, hour, minutes, second), tz


# TODO this should be in the Point class
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


def gaussian(distance, bandwidth):
	"""Calculate a probability that a point originated from a location,
		given the distance."""
	# https://en.wikipedia.org/wiki/Gaussian_function
	a = 1  # height
	b = 0  # mean
	c = bandwidth
	x = distance
	return a * math.exp(-((x-b)**2 / (2*c**2)))


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
	lat1 = math.radians(point1.latitude)
	lat2 = math.radians(point2.latitude)
	diffLong = math.radians(point1.longitude - point2.longitude)
	x = math.sin(diffLong) * math.cos(lat1)
	y = math.cos(lat2) * math.sin(lat1) - (math.sin(lat2) *
		math.cos(lat1) * math.cos(diffLong))
	bearing1 = (math.degrees(math.atan2(x, y))+360) % 360
	# second compass bearing from 2 -> 3
	lat2 = math.radians(point2.latitude)
	lat3 = math.radians(point3.latitude)
	diffLong = math.radians(point3.longitude - point2.longitude)
	x = math.sin(diffLong) * math.cos(lat3)
	y = math.cos(lat2) * math.sin(lat3) -\
		(math.sin(lat2) * math.cos(lat3) * math.cos(diffLong))
	bearing2 = (math.degrees(math.atan2(x, y))+360) % 360
	# we want the smaller of the two angles
	degree_difference = min(abs(bearing1-bearing2),
				(360 - abs(bearing1-bearing2)))
	assert degree_difference <= 180
	return degree_difference


# TODO document and fully understand this algorithm
def viterbi(states, emission_probs, start_probs, transition_probs):
	"""
		states
		emission_probs
		start_probs
		transition_probs
	"""
	V = [{}]
	path = {}
	for state in states:
		# Initialize base cases (t == 0)
		V[0][state] = start_probs[state] * emission_probs[0][state]
		path[state] = [state]
	for t in range(1, len(emission_probs)):	 # Run Viterbi for t > 0
		V.append({})
		newpath = {}
		for s1 in states:
			(prob, state) = max([(V[t-1][s0] * transition_probs[s0][s1] *
				emission_probs[t][s1], s0) for s0 in states])
			V[t][s1] = prob
			newpath[s1] = path[state] + [s1]
		path = newpath  # Don't need to remember the old paths
	(prob, final_state) = max([(V[len(emission_probs)-1][s], s) for s in states])
	# get the optimal sequence of states
	return path[final_state]


def emission_probabilities(point, locations):
	"""
		Given a list of locations, provide the probability that the observed point
		was emitted from none or any. Based on distance and essentially gaussian
		function. The first value in the returned list is the prob that the point
		came from some location not in the list. Right now we're using a guassian
		function that starts 100 meters out from the location.
	"""
	# with d < 100 as at location
	dists = [distance(loc, point) - 100 for loc in locations]
	dists = [0 if d < 0 else d for d in dists]
	probs = [gaussian(d, 100) for d in dists]
	# standardize to one if necessary
	if sum(probs) > 1:
		probs = [p / sum(probs) for p in probs]
	# prepend travel probability as the difference from 1 if there is any
	return [1 - sum(probs)] + probs


def state_transition_matrix(states=[]):
	"""
		Given a list of potential activity location id's, return a simple
		transition probability matrix for use in the viterbi function.
		Transition probs are currently hardcoded and returned as a list of lists.
		0 is the 'travel' state. E.g.:
		____0   1   2  ...
		0  .8  .1  .1
		1  .2  .8  .0
		2  .2  .0  .8
		...
	"""
	# define possible transition probabilities
	travel_to_travel_prob = 0.8
	travel_to_place_prob = 0.2 / len(states)
	place_to_travel_prob = 0.2
	place_to_itself_prob = 0.8
	teleport_prob = 0.0  # impossible, in theory at least
	# make sure travel is an option
	if 0 not in states:
		states += [0]
	trans_prob_matrix = []
	for s0 in states:
		assert type(s0) == int
		trans_prob_matrix.append([])
		for s1 in states:
			if s0 + s1 == 0:  # travel -> travel
				trans_prob_matrix[s0].append(travel_to_travel_prob)
			elif s0 == 0:  # travel -> place
				trans_prob_matrix[s0].append(travel_to_place_prob)
			elif s1 == 0:  # place -> travel
				trans_prob_matrix[s0].append(place_to_travel_prob)
			elif s0 == s1:  # place -> same place
				trans_prob_matrix[s0].append(place_to_itself_prob)
			else:  # place -> place (no travel)
				trans_prob_matrix[s0].append(teleport_prob)
		# TODO I'd like to assert that the rows sum to one, but there 
		# is an annoying amount of floating point error
	return trans_prob_matrix


def read_headers(fname):
	""" (str) -> dict
	Return a dictionary mapping header names to column indices.

	Removes the need to hard coding column numbers when reading files.
	"""
	fd = open(fname)
	d = {}
	header = fd.readline()
	titles = header.split(",")
	for i in range(len(titles)):
		d[titles[i].strip()] = i
	fd.close()
	return d
