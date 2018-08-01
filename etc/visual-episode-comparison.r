USER = 'A'

d = read.csv('~/itinerum/itinerum-trip-breaker/outputs/episodes.csv')
g = read.csv('~/itinerum/itinerum-trip-breaker/outputs/episodes_ground_truth.csv')
# explicitly break into classes: unknown, travel, stationary
d$class = 'travel'
g$class = 'travel'
d[!is.na(d$location_id),'class'] = 'stationary'
g[!is.na(g$location_id),'class'] = 'stationary'
d[d$unknown=='True','class'] = 'unknown'
g[g$unknown=='True','class'] = 'unknown'
d$class = factor(d$class)
g$class = factor(g$class)
# subset to me
d = d[d$user_id == USER,]
g = g[g$user_id == USER,]

# get start and end times
min_start_time = min( c( d$unix_start_time, g$unix_start_time ) )
d$st = (d$unix_start_time - min_start_time) / 3600
g$st = (g$unix_start_time - min_start_time) / 3600
d$et = d$st + c( diff(d$st), 0.1 )
g$et = g$st + c( diff(g$st), 0.1 )

cols = list(stationary='red',travel='green',unknown='blue')

plot( 0, type='n', main=USER, xlim=range(d$st), ylim=c(0,1), bty='n' )
for( i in 1:nrow(d)){ # discovered is on top
	lines(
		x=c( d[i,'st'], d[i,'et'] ),
		y=c( .55, 1 ),
		col=cols[[d[i,'class']]]
	)
}
for( i in 1:nrow(g)){ # ground truth on bottom 
	lines(
		x=c( g[i,'st'], g[i,'et'] ),
		y=c( 0, .45 ),
		col=cols[[g[i,'class']]]
	)
}
