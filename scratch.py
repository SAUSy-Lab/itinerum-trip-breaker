def min_peak(GPS_error_sd,kernel_sd,total_time,threshold_time)
	"""Estimate minimum peak height given time threshold and variance parameters.
		We assume that the volume under the total KDE PDF = 1"""
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
