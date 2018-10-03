# standard modules
import sys, csv
from multiprocessing import Pool, Lock
# our own classes
import config
from points import Location
from trace import Trace


def init_pool():
	global LOCKS
	LOCKS = (Lock(), Lock(), Lock(), Lock())


def initialize_output_files():
	"""Open files for accepting output through script execution."""
	# episodes file
	f = open(config.output_dir+'/episodes.csv', "w+")
	f.write(
		'user_id,sequence,location_id,mode,unknown,local_start_time,'
		'unix_start_time,duration\n')
	f.close()
	# locations file
	f = open(config.output_dir+'/locations.csv', "w+")
	f.write('user_id,location_id,lon,lat,description,used,duration\n')
	f.close()
	# points file
	f = open(config.output_dir+'/classified_points.csv', "w+")
	f.write(
		'user_id,unix_time,known_subset,lon,lat,x,y,weight,removed,timestamp,'
		'interpolated,state,kde,travel_emit_prob\n')
	f.close()
	# days file
	f = open(config.output_dir+'/days.csv', "w+")
	f.write(
		'user_id,date,DoW,'
		'total_time,total_count,'
		'unknown_time,unknown_count,'
		'travel_time,trip_count,'
		'home_time,home_count,'
		'work_time,work_count,'
		'school_time,school_count,'
		'other_time,other_count\n')
	f.close()


def analyze_user(user_data_list):
	"""User data is passed as a list for compatibility with multiprocessing."""
	locks = (None, None, None, None)
	if config.multi_process:
		locks = LOCKS
	user_id = user_data_list[0]
	data = user_data_list[1]
	survey = user_data_list[2]
	user = Trace(user_id, data, survey, locks)
	if len(user.points) > 100:
		if config.debug_output:
			print("User", user_id, 'starts with', len(user.points), 'coordinates')
		user.remove_repeated_points()
		user.remove_known_error(config.min_accuracy)
		user.remove_sequential_duplicates()
		user.remove_positional_error()
		# this is actually necessary again after positional cleaning
		# ( some angles == 0 )
		user.remove_sequential_duplicates()
		# identify gaps in the data
		user.make_known_subsets()
		if len(user.known_subsets) > 0:
			# find locations with the cleaned data
			user.get_activity_locations()
			# allocate time
			user.break_trips()
			user.identify_named_locations()
			# write the output
			user.flush()
		else:  # not enough subsets of appropriate length
			print("\tinsufficient data for {}".format(user_id))


if __name__ == "__main__":
	# Get all the coordinates data in a big list per user so that we only
	# have to read this file once.
	user_data = {}
	# get a list of all users in the coordinates file
	with open(config.input_coordinates_file, newline='') as f:
		reader = csv.DictReader(f)
		for row in reader:
			if row['uuid'] not in user_data:
				user_data[row['uuid']] = [row]
			else:
				user_data[row['uuid']].append(row)
	# Check for the existence of user provided location data from the survey 
	# read in once up front and pass as a dict
	user_locations = {}
	with open(config.input_survey_responses_file, newline='') as f:
		reader = csv.DictReader(f)
		for record in reader:
			user_id = record['uuid']
			user_locations[user_id] = {}
			if 'location_home_lon' in record and record['location_home_lon'] != '':
				user_locations[user_id]['home'] = Location(
					record['location_home_lon'], record['location_home_lat'] 
				)
			if 'location_work_lon' in record and record['location_work_lon'] != '':
				user_locations[user_id]['work'] = Location(
					record['location_work_lon'], record['location_work_lat']
				)
			if 'location_study_lon' in record and record['location_study_lon'] != '':
				user_locations[user_id]['school'] = Location(
				record['location_study_lon'], record['location_study_lat']
			)
	if config.debug_output:
		print(len(user_data), 'user(s) to clean')
	# loop over users calling all the functions for each
	initialize_output_files()
	# create a list of users' data to work on
	list_of_users = []
	for uid, data in user_data.items():
		# uncomment the following line for testing a single user
		if uid != '': continue
		list_of_users.append((uid, data, user_locations[uid]))
	# parallel processing option
	if config.multi_process:
		LOCKS = (Lock(), Lock(), Lock(), Lock())
		p = Pool(processes=config.num_pro, initializer=init_pool)
		p.map(analyze_user, list_of_users)
	# non-parallel processing
	else:
		for user_data in list_of_users:
			analyze_user(user_data)

	print("Done!", file=sys.stderr)
