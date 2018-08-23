# functions related to Guassian kernel functions or normal distributions

from math import exp
import config
from scipy.stats import multivariate_normal


def gaussian(x, bandwidth):
	"""
	Evaluate a Gaussian function at x for a distribution centered on zero
	with a height of 1 given a bandwidth.
	Used as a distance decay function.
	See: https://en.wikipedia.org/wiki/Gaussian_function
	"""
	return exp(-(x**2 / (2 * bandwidth**2)))


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
	if config.debug_output:
		print('\tRunning KDE on', len(x_vector), 'points')
	point_matrix = cbind(FloatVector(x_vector), FloatVector(y_vector))
	bandwidth = config.kernel_bandwidth
	surface = ks.kde(x=point_matrix, eval_points=point_matrix,
			w=FloatVector(weights), H=diag(FloatVector([bandwidth**2, bandwidth**2])),
binned=False)
	estimates = surface.rx2('estimate')
	# turn these into more pythonish objects so that the rpy2 syntax doesn't
	# have to leave this function
	est = []
	for i in range(1, len(weights)+1):
		# insert estimate values
		est.append(estimates.rx(i)[0])
	# this is now a vector (python list)
	return est
