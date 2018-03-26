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


def clean_sequence(sequence):
	"""DOCUMENTATION NEEDED"""
	pass


def init_file(filename, t):
	"""DOCUMENTATION NEEDED"""
	fd = open(filename, "w")
	header = ""
	if t == "activities":
		header = "user_id,sequence,location_id,travel_mode(s),Unknown,start_time\n"
	elif t == "locations":
		header = "user_id,uid,lon,lat,description\n"
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
	init_file(config.FILENAME, "activities")
	for user_id in user_ids:
		user = Trace(user_id)
		print( len(user.points),'points at start for',user_id )
		user.remove_known_error(100)
		user.remove_sequential_duplicates()
		user.remove_positional_error()
		# this is actually necessary again after positional cleaning
		# ( some angles == 0 )
		user.remove_sequential_duplicates()

		# GEOTESTING: checking post-cleaning geometry
		import csv
		with open('outputs/TESTING_post-cleaning-points.csv', 'w+') as csvfile:
			writer = csv.writer(csvfile, delimiter=',', quotechar='"')
			writer.writerow(['latitude','longitude','timestamp'])
			for point in user.points:
				writer.writerow([point.latitude,point.longitude,point.time])

		user.make_subsets()
		user.break_trips() 


