def min_peak(GPS_error_sd,kernel_sd,total_time,threshold_time):
	"""Estimate minimum peak height given time threshold and variance parameters.
		We assume that the volume under the total KDE PDF ~= 1. 
		-	total_time is the sum of the time weights used in the KDE.
		-	threshold_time is the minimum activity time."""
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
	return peak_height * (threshold_time / total_time)


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
	"""Inputs are 2D spatial grids where cells are indexed by consecutive 
		integers e.g. grid[column][row].
		Estimates gives the estimated probability at a point
		Locations gives a geographic location for that estimate
		Threshold is a minimum height that a peak must reach.
		Return a list of lon,lat tuples for potential activity locations"""
	print( '\tFinding peaks in PDF surface' )
	def get_starting_cell():
		"""get the indices of a cell with a value of one"""
		for x, column in enumerate(to_visit):
			for y, value in enumerate(column):
				if value == 1:
					return (x,y)
		assert False # should never be here
	def span_cluster(x,y):
		"""Breadth-first search function. Visit all contiguous cells with a 
			value of one (from an initial point) and return a list of their 
			x,y indices."""
		queue = [(x,y)]
		discovered = []
		# set limits for later
		max_x = len(to_visit)-1
		max_y = len(to_visit[0])-1
		# while there is stuff in the queue
		while len(queue) > 0:
			# visit the first thing in
			x,y = queue.pop(0)
			# mark it discovered
			discovered.append((x,y))
			# check all possible neighbor locations (rook contiguity)
			for delta_x,delta_y in [(0,-1),(-1,0),(+1,0),(0,+1)]:
				new_x = x + delta_x
				new_y = y + delta_y
				# check for any reason not to visit this cell
				if (new_x,new_y) in discovered: continue 
				if (new_x,new_y) in queue: continue
				if new_x < 0 or new_y < 0: continue
				if new_x > max_x or new_y > max_y: continue
				if to_visit[new_x][new_y] != 1: continue
				# everything checked out and this cell is part of the cluster
				queue.append((new_x,new_y))
		# update the matrix that this cluster has all beeen visited
		for x,y in discovered:
			to_visit[x][y] = 0
		return discovered
	# Create an array of the same dimensions to keep track of cells we need to 
	# visit. 1 for a yet unvisited cells, 0 otherwise. We make the simplifying 
	# assumption that we don't need to visit cells under the threshold.
	to_visit = [ [ 1 if value > threshold else 0 for value in column ] for column in estimates ]
	# now detect clusters of cells
	# while there are still true cells unvisited
	clusters = []
	# cells are set = to 0 once visited
	while sum([sum(row) for row in to_visit]) > 0:
		# get an arbitrary starting point
		x,y = get_starting_cell()
		print(x,y, sum([sum(row) for row in to_visit]))
		# call a recursive function to visit all neighbors
		clusters.append( span_cluster(x,y) )
	print( '\t',len(clusters),'clusters found' )
	peaks = []
	# each 'cluster' is a list of cells in a cluster
	for cluster in clusters:
		# find the peak as the maximum value of cells in the cluster
		values = [ estimates[column][row] for column,row in cluster ]
		cluster_max = max(values)
		# now find the location of the peak
		for column,row in cluster:
			if estimates[column][row] == cluster_max:
				peaks.append(locations[column][row])
				break

	# GEOTESTING: checking sampling geometry
	import csv
	with open('outputs/TESTING_potential-activity-locations.csv', 'w+') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quotechar='"')
		writer.writerow(['longitude','latitude'])
		for x,y in peaks:
			lon,lat = unproject(x,y))
			writer.writerow([lon,lat])




def kde(x_vector,y_vector,weights,bandwidth,cell_size):
	"""Do weighted 2d KDE in R KS package, returning python results.
		Returns two 2d arrays: estimates and estimate locations."""
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
	from rpy2.robjects import FloatVector, IntVector
	# set the range to the bounding box plus some
	x_min = min(x_vector) - bandwidth
	x_max = max(x_vector) + bandwidth
	y_min = min(y_vector) - bandwidth
	y_max = max(y_vector) + bandwidth
	# set the approximate cell size by determining the number of cells from it
	num_cells_x = int( (x_max-x_min)/cell_size )
	num_cells_y = int( (y_max-y_min)/cell_size )
	print( '\tRunning KDE on',num_cells_x,'x',num_cells_y,'grid with',len(x_vector),'points' )
	# do the KDE
	surface = ks.kde(
		x = cbind( FloatVector(x_vector), FloatVector(y_vector) ),
		H = diag( FloatVector( [ bandwidth**2, bandwidth**2 ] ) ),
		xmin = FloatVector( [x_min,y_min] ),
		xmax = FloatVector( [x_max,y_max] ),
		gridsize = IntVector( [num_cells_x,num_cells_y] )
	)
	eval_points = surface.rx2('eval.points')
	estimates = surface.rx2('estimate')
	# turn these into more pythonish objects so that the rpy2 syntax doesn't 
	# have to leave this function
	eva, est = [], []
	for x in range(1,num_cells_x+1):
		est.append([])
		eva.append([])
		for y in range(1,num_cells_y+1):
			# insert estimate values
			est[x-1].append(estimates.rx(x,y)[0])
			# insert location tuples
			easting, northing = eval_points.rx(1)[0][x-1], eval_points.rx(2)[0][y-1]
			eva[x-1].append( (easting, northing) )
	# these are now 2d arrays (python lists) giving estimated probabilities
	# and locations as x,y tuples
	# both lists are indexed as list[x][y]
	assert len(est) == num_cells_x
	assert len(est[1]) == num_cells_y

	# GEOTESTING: checking sampling geometry
	import csv
	with open('outputs/TESTING_kde-eval-points.csv', 'w+') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quotechar='"')
		writer.writerow(['x','y','estimate'])
		for x, column in enumerate(eva):
			for y, (east,north) in enumerate(column):
				writer.writerow([east,north,est[x][y]])

	return est, eva

