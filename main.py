# This is an early algorithm for cleaning itinerum output data.
# I'm mainly trying to remove obvious positional errors and drop out 
# redundant points. After this comes trip/activity partitioning. 
# 
# The key thing here is to remove points that don't appear to be based 
# on reasonably accurate *GPS*. Id est, some points come from black-box 
# systems inside the phone, and some are GPS based but not accurate. 
# There a few things we can look for to pick these out:

# 1) GPS noise has a roughly gaussian distribution, and should very 
#		rarely hit the same exact point twice, given enough precision. 
#		Yet phones often get "stuck" on a point and repeat it. 
#		This may indicate a cell-tower or some other problem. Any 
#		point that is used two or more times should come under strict 
#		scrutiny. This is not quite yet implemented here. 

# 2) Any major jump away and back again, especially if the h_error
#		Doesn't justify the distance may indicate a non-GPS signal 
#		or a bad error estimate. Away-and-back-again points are identified 
#		by the minimum distance from neighboring points and the angle formed 
#		between the three. 

# The algorithm iterates over users and advances iteratively over each 
# user's data.  

# standard modules
import datetime, csv, math, rpy2
# our own classes
from point import Point
from trace import Trace
import config


def init_file(filename, t):
	"""DOCUMENTATION NEEDED"""
	fd = open(filename, "w")
	header = ""
	if t == "activities":
		header = "user_id,sequence,location_id,travel_mode(s),Unknown,start_time\n"
	elif t == "locations":
		header = "user_id,uid,lon,lat,description,time_at\n"
	fd.write(header)
	fd.close()


# Standard format so we can import this module elsewhere.
if __name__ == "__main__":
	user_ids = []
	# get a list of all user_ids in the coordinates file
	with open(config.input_coordinates_file, newline='') as f:
		reader = csv.DictReader(f)
		for row in reader:
			user_ids.append( row['uuid'] )
	# Keep only unique user_id's
	user_ids = list(set(user_ids)) 
	print( len(user_ids),'user(s) to clean' )
	# loop over users calling all the functions for each
	init_file(config.output_activities_file, "activities")
	init_file(config.output_locations_file, "locations")
	u = 1
	for user_id in user_ids:
		# create trace object for this user
		user = Trace(user_id)
		# remove GPS points believed to be in error
		print("User :", u, len(user.points),'points at start for',user_id )
		u += 1
		user.remove_known_error( config.min_accuracy )
		user.remove_sequential_duplicates()
		user.remove_positional_error()
		# this is actually necessary again after positional cleaning
		# ( some angles == 0 )
		user.remove_sequential_duplicates()
		# identify gaps in the data
		user.make_known_subsets()
		# find locations with the cleaned data
		user.get_activity_locations()
		# allocate time
		user.break_trips() 


