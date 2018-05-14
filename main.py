# standard modules
import datetime, csv, math, rpy2
# our own classes
from point import Point
from trace import Trace
import config

def initialize_output_files():
	"""Open files for accepting output through script execution."""
	# episodes file
	f = open(config.output_episodes_file, "w")
	f.write('user_id,sequence,location_id,mode,unknown,start_time\n')
	f.close()
	# locations file
	f = open(config.output_locations_file, "w")
	f.write('user_id,location_id,lon,lat,description,used\n')
	f.close()
	# points file
	f = open(config.output_points_file, "w")
	f.write('user_id,lon,lat,weight,removed,interpolated,state\n')
	f.close()
	# days file
	f = open(config.output_days_file, "w")
	# TODO so much more to do here
	f.write('user_id,date,DoW,total_minutes,trip_count,travel_time,unknown_time,home_time,work_time,school_time,home_count,work_count,school_count\n')
	f.close()

# Standard format so we can import this module elsewhere.
if __name__ == "__main__":
	user_ids = {}
	# get a list of all user_ids in the coordinates file
	with open(config.input_coordinates_file, newline='') as f:
		reader = csv.DictReader(f)
		for row in reader:
			if row['uuid'] not in user_ids:
                        	user_ids[row['uuid'] = [row]
			else:
				user_ids[row['uuid'].append(row)
	print( len(user_ids),'user(s) to clean' )
	# loop over users calling all the functions for each
	initialize_output_files()
	u = 1
	for user_id in user_ids:
		try:
			# create trace object for this user
			user = Trace(user_id, user_ids[user_ids])
			if len(user.points) < 100: continue
			# remove GPS points believed to be in error
			print("User:", u, len(user.points),'points at start for',user_id )
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
			# 
			user.identify_locations()
			# allocate time
			user.break_trips()
			# write the output
			user.flush()
		except Exception: # TODO specify what exceptions are permissible
			print( user.id,'aborted' )


