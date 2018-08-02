USER = 'A'

cols = c(
		rgb(1,0,0,alpha=0.5), # red: stationary
		rgb(0,1,0,alpha=0.5), # green: travel
		rgb(0,0,1,alpha=0.5)  # blue: unknown
	)

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
d$st = d$unix_start_time
g$st = g$unix_start_time
d$et = d$unix_start_time + c( diff(d$st), 0.1 )
g$et = g$unix_start_time + c( diff(g$st), 0.1 )

pdf(paste0('~/',USER,'.pdf'),width=3,height=200)
	par(
		# bottom left top right
		mar=c(0,0,0,0),
		mai=c(0,1.4,0,0),
		family='mono'
	) 
	min_hour = min(d$st) - min(d$st) %% 3600
	max_hour = max(d$st) + max(d$st) %% 3600
	plot( 0, type='n', 
		main=USER, 
		ylim=c(min_hour,max_hour), xlim=c(0,1), 
		bty='n', xaxt='n', yaxt='n',xlab='',ylab=''
	)
	axis( 2, at=seq(min_hour,max_hour,3600), las=2 )
	rect( # discovered is on right
		xleft=0.55, xright=1,
		ytop=d$st, ybottom=d$et,
		col=c(cols[d$class])
	)
	rect( # ground truth is on left
		xleft=0, xright=.45,
		ytop=g$st, ybottom=g$et,
		col=c(cols[g$class])
	)
dev.off()