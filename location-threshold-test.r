# script to determine the height of a density surface for a stationary activity
# given variable GPS error
# WORK IN PROGRESS

# stdev of distance of point from actual location
GPS_error = 10
# cell size in meters
cell_size = 0.5
# num points
n = 1000
# bandwidth in meters
bandwidth = 10

xmin = c(-100,-100)
xmax = c(100,100)
num_cells_x = abs( xmin[1] - xmax[1] ) / cell_size
num_cells_y = abs( xmin[2] - xmax[2] ) / cell_size
x = rnorm( n, sd=GPS_error )
y = rnorm( n, sd=GPS_error )

library(ks)

surface = kde(
	x = cbind(x,y), 
	H=matrix(c(bandwidth,0,0,bandwidth),2),
	gridsize = c(num_cells_x,num_cells_y),
	xmin = xmin,
	xmax = xmax,
	verbose = TRUE
)
# whats the maximum of the estimated probabilities
print( 'est stdev of joint distribution' )
new_sd = sqrt( GPS_error^2 + bandwidth^2 )
print(new_sd)
print( 'observed peak of surface' )
print(max(surface$estimate))
print( 'theoretical peak of surface' )
# estimate the bivariate normal PDF at 0,0 
est_peak = dmvnorm( x=c(0,0), mean=c(0,0), sigma=diag(c(new_sd,new_sd)) )
print( est_peak  )
