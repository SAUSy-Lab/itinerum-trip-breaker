# functions related to Guassian kernel functions or normal distributions

from math import exp
import config
from scipy.stats import multivariate_normal


def gaussian_log_prob(x, bandwidth):
	""" Evaluate a Gaussian function at x for a distribution centered on zero
	with a height of 1 given a bandwidth. Return the natural log of the 
	probability. Used as a distance decay function.
	See: https://en.wikipedia.org/wiki/Gaussian_function
	"""
	return -(x**2 / (2 * bandwidth**2))
	# unlogged: return exp(-(x**2 / (2 * bandwidth**2)))


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


def kde(input_points,eval_points):
	"""Do weighted 2d KDE in R KS package, returning python results.
	Returns list of P estimates corresponding to estimate_points.
	input_points should be a list of points with x,y,weight attributes
	eval_points are points where the PDF will be estimated (not elsewhere).	

	Another possible way of doing KDE (in Python) is with
	http://pysal.readthedocs.io/en/latest/users/tutorials/
	smoothing.html#non-parametric-smoothing ?
	or with http://scikit-learn.org/stable/modules/density.html ?

	But I'm not convinced yet that it's worthwhile to explore these.
	R is quite fast and works fine. 
	"""
	# import R tools
	from rpy2.robjects.packages import importr  # for importing packages
	ks = importr('ks') # load the ks package (making kde function available)
	from rpy2.robjects import r, FloatVector  # variable names
	cbind, diag = r['cbind'], r['diag'] # get some basic R functions into Python
	# create the inputs from input_points
	x_in, y_in = [p.x for p in input_points], [p.y for p in input_points]
	weights  = [p.weight for p in input_points]
	# standardize the weights to the sample size
	if sum(weights) != len(weights):
		adjust_factor = len(weights) / float(sum(weights))
		weights = [w * adjust_factor for w in weights]
	# get eval_points processed for R
	x_ev, y_ev = [p.x for p in eval_points], [p.y for p in eval_points]
	# do the KDE
	if config.debug_output:
		print('\tRunning KDE on', len(input_points), 'points')
	bandwidth = config.kernel_bandwidth
	surface = ks.kde(
		x = cbind(FloatVector(x_in), FloatVector(y_in)),
		eval_points = cbind(FloatVector(x_ev), FloatVector(y_ev)),
		w = FloatVector(weights),
		H = diag( FloatVector([bandwidth**2, bandwidth**2]) ),
		binned = False
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
