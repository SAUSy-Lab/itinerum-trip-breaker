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

def find_peaks(grid,threshold):
	"""this is sort of a placeholder function which will need to be written in 
		an implementation specific way once we have KDEs calculated. 
		Input here is a 2D spatial grid where cells are indexed by consecutive 
		integers e.g. grid[row][column].
		Threshold is a minimum height that a peak must reach."""
	# create a copy with binary threshold values
	from copy import copy
	bin_grid = copy(grid)
	# compare to threshold
	for row_id, column in enumerate( bin_grid ):
		for column_id, cell_value in enumerate(bin_grid[row_id]):
			cell_value = 1 if cell_value >= threshold else 0
	# now detect clusters of True cells
	# while there are still true cells unvisited
	clusters = []
	while sum(bin_grid) > 0:
		# recursively visit neighboring cells and add their locations to a list 
		# if they are past the threshold. Once visited, set value to zero.
		pass
	peaks = []
	# each 'cluster' is a list of cells in a cluster
	for cluster in clusters:
		# find the peak as the maximum value of cells in the cluster
		values = [ grid[row][column] for row,column in cluster ]
		cluster_max = max(values)
		# now find the lcation of the peak
		for row, column in cluster:
			if grid[row][column] == cluster_max:
				peaks.append(row,column)
				break
		# peaks should now be the cell locations (row,column) of the maxima of 
		# each cluster. These will need to be mapped back to geographical space.
		# we will know the cell size and the origin of the grid, so this should 
		# be pretty easy. These points are our activity locations.


def kde(x_vector,y_vector,weights,bandwidth,cell_size=1):
	"""Do weighted 2d KDE in R KS package, returning python results.
		Returns two 2d arrays: estimates and estimate locations."""
	# check the inputs
	assert len(x_vector) == len(y_vector)
	assert len(weights) == len(x_vector)
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
	x_min = min(x_vector) - bandwidth * 2
	x_max = max(x_vector) + bandwidth * 2
	y_min = min(y_vector) - bandwidth * 2
	y_max = max(y_vector) + bandwidth * 2
	# set the cell size to roughly 1 meter 
	num_cells_x = int( (x_max-x_min)/cell_size )
	num_cells_y = int( (y_max-y_min)/cell_size )
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
	assert len(est) == len(num_cells_x)
	assert len(est[1]) == len(num_cells_y)
	return est, eva
	

## testing
#x = [1,2,3,4.5,5]
#y = [1,2,3.5,4,6]
#w = [1,2,1,2,2]
#b = 5
#estimates, eval_points = kde(x,y,w,b)

#print 'top left cell estimate', estimates.rx(1,1)[0]
#print 'top left cell location', eval_points.rx(1)[0][0],',',eval_points.rx(2)[0][0]
#print 'row_length',len(estimates.rx(1,True))
##print eval_points
##print sum(estimates)

