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


def find_peaks(estimates,locations,threshold):
	"""PDF was estimated at a selection of points, which are here given as a list
		of P values (estimates) and a list of (x,y) locations. The idea is to toss 
		out any values below the threshold and identify spatial clusters among 
		those that remain. In each such cluster, the highest value is the activity 
		location."""
	assert len(estimates) == len(locations)
	# drop values below the threshold
	locations = [ (x,y) for (x,y),est in zip(locations,estimates) if est >= threshold ]
	estimates = [ est for est in estimates if est >= threshold ]
	assert len(estimates) == len(locations)
	print('\tclustering',len(estimates),'points above',threshold,'threshold')


#def find_peaks_breadth_first(estimates,locations,threshold):
#	"""Inputs are 2D spatial grids where cells are indexed by consecutive 
#		integers e.g. grid[column][row].
#		Estimates gives the estimated probability at a point
#		Locations gives a geographic location for that estimate
#		Threshold is a minimum height that a peak must reach.
#		Return a list of lon,lat tuples for potential activity locations"""
#	print( '\tFinding peaks in PDF surface' )
#	def get_starting_cell():
#		"""get the indices of a cell with a value of one"""
#		for x, column in enumerate(to_visit):
#			for y, value in enumerate(column):
#				if value == 1:
#					return (x,y)
#		assert False # should never be here
#	def span_cluster(x,y):
#		"""Breadth-first search function. Visit all contiguous cells with a 
#			value of one (from an initial point) and return a list of their 
#			x,y indices."""
#		queue = [(x,y)]
#		discovered = []
#		# set limits for later
#		max_x = len(to_visit)-1
#		max_y = len(to_visit[0])-1
#		# while there is stuff in the queue
#		while len(queue) > 0:
#			# visit the first thing in
#			x,y = queue.pop(0)
#			# mark it discovered
#			discovered.append((x,y))
#			# check all possible neighbor locations (rook contiguity)
#			for delta_x,delta_y in [(0,-1),(-1,0),(+1,0),(0,+1)]:
#				new_x = x + delta_x
#				new_y = y + delta_y
#				# check for any reason not to visit this cell
#				if (new_x,new_y) in discovered: continue 
#				if (new_x,new_y) in queue: continue
#				if new_x < 0 or new_y < 0: continue
#				if new_x > max_x or new_y > max_y: continue
#				if to_visit[new_x][new_y] != 1: continue
#				# everything checked out and this cell is part of the cluster
#				queue.append((new_x,new_y))
#		# update the matrix that this cluster has all beeen visited
#		for x,y in discovered:
#			to_visit[x][y] = 0
#		return discovered
#	# Create an array of the same dimensions to keep track of cells we need to 
#	# visit. 1 for a yet unvisited cells, 0 otherwise. We make the simplifying 
#	# assumption that we don't need to visit cells under the threshold.
#	to_visit = [ [ 1 if value > threshold else 0 for value in column ] for column in estimates ]
#	# now detect clusters of cells
#	# while there are still true cells unvisited
#	clusters = []
#	# cells are set = to 0 once visited
#	while sum([sum(row) for row in to_visit]) > 0:
#		# get an arbitrary starting point
#		x,y = get_starting_cell()
#		print(x,y, sum([sum(row) for row in to_visit]))
#		# call a recursive function to visit all neighbors
#		clusters.append( span_cluster(x,y) )
#	print( '\t',len(clusters),'clusters found' )
#	peaks = []
#	# each 'cluster' is a list of cells in a cluster
#	for cluster in clusters:
#		# find the peak as the maximum value of cells in the cluster
#		values = [ estimates[column][row] for column,row in cluster ]
#		cluster_max = max(values)
#		# now find the location of the peak
#		for column,row in cluster:
#			if estimates[column][row] == cluster_max:
#				peaks.append(locations[column][row])
#				break
#	# GEOTESTING: checking sampling geometry
#	import csv
#	with open('outputs/TESTING_potential-activity-locations.csv', 'w+') as csvfile:
#		writer = csv.writer(csvfile, delimiter=',', quotechar='"')
#		writer.writerow(['longitude','latitude'])
#		for x,y in peaks:
#			lon,lat = unproject(x,y)
#			writer.writerow([lon,lat])




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

