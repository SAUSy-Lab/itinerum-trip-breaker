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
	CLUSTER_DISTANCE = 50 # meters
	from math import sqrt
	# drop values below the threshold
	locations = [ (x,y) for (x,y),est in zip(locations,estimates) if est >= threshold ]
	estimates = [ est for est in estimates if est >= threshold ]
	assert len(estimates) == len(locations)
	print('\tclustering',len(estimates),'points above',threshold,'threshold')
	# now calculate a distance-based connectivity matrix between all these points
	neighbs = []
	for i,(x1,y1) in enumerate(locations):
		neighbs.append([])
		for j,(x2,y2) in enumerate(locations):
			# use euclidian distance since this is already projected
			connection = sqrt((x1-x2)**2 + (y1-y2)**2) < CLUSTER_DISTANCE
			neighbs[i].append( connection )
	print( '\thave connection matrix with',sum([len(l) for l in neighbs]),'entries' )
	# clusters will be a list of disjoint sets
	clusters = []
	# now for each point, check for cluster membership and add and merge clusters
	for i in range(0,len(neighbs)):
		# get a set of neighbor indices within distance, 
		# including this point itself ( dist = 0 )
		neighbors = set( [ j for j,n in enumerate(neighbs[i]) if n ] )
		# create list to keep track of clusters this set belongs to
		member_cluster_indices = []
		for i,cluster in enumerate(clusters):
			# check each cluster for overlap with this set
			if not cluster.isdisjoint( neighbors ):
				member_cluster_indices.append(i)
		if len(member_cluster_indices) == 0:
			#  we have no overlap, so this becomes a new cluster
			clusters.append(neighbors)
		elif len(member_cluster_indices) > 0:
			# we have one or more matching clusters
			for i in reversed(member_cluster_indices):
				# add everyhting together in a new cluster and
				# drop off the old clusters which are now merged in the new one
				neighbors = neighbors | clusters.pop(i)
				# add the new cluster
			clusters.append(neighbors)

	print( '\tfound',len(clusters),'clusters with',sum([len(c) for c in clusters]),'total points' )
	from location import ActivityLocation
	potential_activity_locations = []
	# find the maximum estimate and a location with that value
	for cluster in clusters:
		peak_height = max( [estimates[i] for i in cluster] )
		for i in cluster:
			# if this is the highest point
			if estimates[i] == peak_height:
				x,y = locations[i]
				lon,lat = unproject(x,y)
				# create a location and append to the list
				location = ActivityLocation(lon,lat)
				potential_activity_locations.append( location )
				break

	# GEOTESTING: checking post-cleaning geometry
	import csv
	with open('outputs/TESTING_potential-activity-locations.csv', 'w+') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quotechar='"')
		writer.writerow(['longitude','latitude'])
		for location in potential_activity_locations:
			writer.writerow([location.longitude,location.latitude])

	return potential_activity_locations



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

