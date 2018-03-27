# script to verify the height of a density surface for a stationary activity
# given variable GPS error and kernel bandwidth.
#
# The idea behind this script is to verify that we can analytically determine the 
# height of a peak in the estimated PDF given a known GPS error and kernel bandwidth.
# I consider this so verified. - Nate

# stdev of distance of point from actual location
GPS_error = 50
# bandwidth in meters (standard_deviation)
bandwidth = 50
# cell size in meters
cell_size = 0.5
# num points
n = 1000

xmin = c(-100,-100)
xmax = c(100,100)
num_cells_x = abs( xmin[1] - xmax[1] ) / cell_size
num_cells_y = abs( xmin[2] - xmax[2] ) / cell_size
x = rnorm( n, sd=GPS_error )
y = rnorm( n, sd=GPS_error )

library(ks)

surface = kde(
	x = cbind(x,y), 
	H=diag(c(bandwidth^2,bandwidth^2)), # this is a covariance matrix
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
sd_m = diag( c( new_sd^2, new_sd^2 ) )
est_peak = dmvnorm( x=c(0,0), mean=c(0,0), sigma=sd_m )
print( est_peak )

# plot x and y cross sections of the plot
plot(
	y = surface$estimate[200,],
	x = surface$eval.points[[1]],
	type = 'l',
	col='red',
	main='bivariate normal cross sections at 0,0',
	xlab='PDF estimate',
	ylab='x/y'
)
lines(
	y = surface$estimate[,200],
	x = surface$eval.points[[2]],
	type = 'l',
	col='blue'
)
# draw a line where the peak should be based on the analytical solution
abline(h=est_peak)
