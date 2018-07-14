# standard modules
import datetime
import sys
import csv
import math
import rpy2
# our own classes
from point import Point
from trace import Trace
from location import Location
from multiprocessing import Pool
import config


def initialize_output_files():
	"""
	Open files for accepting output through script execution.
	"""
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
	f.write('user_id,lon,lat,x,y,weight,removed,interpolated,state,kde\n')
	f.close()
	# days file
	f = open(config.output_days_file, "w")
	f.write('user_id,date,DoW,total_minutes,trip_count,travel_time,')
	f.write('unknown_time,home_time,work_time,school_time,home_count,')
	f.write('work_count,school_count\n')
	f.close()


def analyze_user(tup):
	user_id = tup[0]
	data = tup[1]
	survey = tup[2]
	user = Trace(user_id, data, survey)
	if len(user.points) > 100:
		if config.db_out:
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
		# find locations with the cleaned data
		user.get_activity_locations()
		user.remove_repeated_locations()
		# allocate time
		user.break_trips()
		# write the output
		user.flush()


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
	# read in the survey response location data once now and store in a dict
	survey_responses = {}
	with open(config.input_survey_responses_file, newline='') as f:
		reader = csv.DictReader(f)
		for row in reader:
			home = Location(row['location_home_lon'], row['location_home_lat'])
			work = Location(row['location_work_lon'], row['location_work_lat'])
			if row['location_study_lon'].strip() and row['location_study_lat'].strip():
				school = Location(row['location_study_lon'].strip(),
				row['location_study_lat'].strip())
			else:
				school = None
			survey_responses[row['uuid']] = [home, work, school]
	if config.db_out:
		print(len(user_data), 'user(s) to clean')
	# loop over users calling all the functions for each
	initialize_output_files()
	# loop over users calling all the functions for each

	list_of_users = []
	for uid, data in user_data.items():
		# survey responses[ui] instead of (None, None, None)
		list_of_users.append((uid, data, (None, None, None)))
	if config.multi_process:
		with Pool(config.num_pro) as p:
			p.map(analyze_user, list_of_users)
	else:
		for t in list_of_users:
			analyze_user(t)
	print("Done!\7", file=sys.stderr)
