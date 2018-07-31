from point import Point
from config import *
from misc_funcs import distance
from statistics import median
from datetime import timedelta, datetime
from misc_funcs import parse_ts, read_headers

# Location comparison functions

def compare_locations(truth, compd):
	""" (str, str) -> [(str, int, float, float)]
	Compare ground truth locations to the computed locations,
	and output a list of quality metrics by user.
	"""
	true_locations = get_locations(truth)
	computed_locations = get_locations(compd)
	
	users_to_matrix = {}
	results = []
	for user in true_locations:
		if user not in computed_locations:
			print(user, " not in computed locations")
		else:
			num_locations = len(computed_locations[user]) - len(true_locations[user])
			num_to_match = min(len(computed_locations[user]), len(true_locations[user]))
			users_to_matrix[user] = {}
			# TODO exclude unused locations
			h = read_headers(truth)
			distance_matrix(h, users_to_matrix[user],
					true_locations[user], computed_locations[user])
			# users_to_matrix[user][true_location][computed_location] = distance
			# dictionary of dictionaries of dictionaries of distances

			min_distances = compare_user_locations(users_to_matrix[user], num_to_match)
			mean_min_dis = sum(min_distances) / len(min_distances)
			med = median(min_distances)
			results.append((user, num_locations, mean_min_dis, med))
	return results

def compare_user_locations(distances, num_to_match):
	min_distances = []
	dist_list = []
	for loc in distances:
		for guess in distances[loc]:
			dist_list.append((distances[loc][guess], guess, loc))
	for _ in range(num_to_match):
		dist_list.sort()
		best = dist_list.pop(0)  # First in the list has the lowest distance
		min_distances.append(best[0])  # save the distance entry
		compd = best[1]  # the computed location
		truth = best[2]  # the true location
		for entry in dist_list:
			# if the computed or true location match
                        # remove the rest of the entries with those locations
			if entry[1] == compd or entry[2] == truth:
				dist_list.remove(entry)
	return min_distances


def distance_matrix(h, matrix, truths, compds):
	""" (dict, dict, dict, dict) -> None
	Update matrix for user with a 2d matrix mapping true
	locations and computed locations to their respective distances.
	"""
	for location in truths:
		matrix[location[h["location_id"]]] = {}
		for guess in compds:
			p1 = Point("", guess[h["lon"]], guess[h["lat"]], 0)
			p2 = Point("", location[h["lon"]], location[h["lat"]], 0)
			# Felipevh forgets if these need to be projected first
			matrix[location[h["location_id"]]][guess[h["location_id"]]] = distance(p1, p2)

def get_locations(location_file):
	""" (str) -> dict
	Return a dictionary mapping users to the location described in 
	location_file
	"""
	utl = {}  # users to locations dictionary
	fd = open(location_file)
	h = read_headers(location_file)
	line = fd.readline().split(',')
	for line in fd:
		loc = line.split(',')
		if loc[h["user_id"]] not in utl:
			utl[loc[h["user_id"]]] = [loc]
		else:
			utl[loc[h["user_id"]]].append(loc)
	return utl

# Episode comparison functions below

def compare_episodes(truth, compd):
	"""
	"""
	true_users_to_eps = read_episodes(truth)
	computed_users_to_eps = read_episodes(compd)
	# get the headers, both files should match
	h = read_headers(truth)
	metrics = []
        
	for user in true_users_to_eps:
		if user in computed_users_to_eps:
			user_metrics = (user,) + compare_user_episodes(true_users_to_eps[user], computed_users_to_eps[user], h)
			metrics.append(user_metrics)
		else:
			print("{} not in computed episodes".format(user))
	return metrics

def read_episodes(file_name):
	"""
	"""
	h = read_headers(file_name)
	fd = open(file_name)
	users_to_eps = {}
	fd.readline()
	for l in fd:
		line = l.strip().split(",")
		user = line[h["user_id"]]
		if user in users_to_eps:
			users_to_eps[user].append(line)
		else:
			users_to_eps[user] = [line]
	fd.close()
	# Sort the lists of episodes by unix time
	for u in users_to_eps:
		# the lambda expression sorts the episode entries by their unix time
		users_to_eps[u].sort(key=lambda x: x[-1])	
	return users_to_eps

def compare_user_episodes(true, computed, headers):
	end_time = min(true[-1][-1], computed[-1][-1])
	start_time = max(true[0][-1], computed[0][-1])
	i, j = 0, 0
	# First we bring indices to overlapping episodes
	while true[i+1] < start_time:
		i = i + 1
	while computed[j+1] < start_time:
		j = j + 1

	# Iterate over episodes incrementally so that they always overlap
	while true[i] < end_time and computed[j] < end_time:
		pass  # do work
	return (1,2,3,4)

def literal_eval(string):
	""" (str) -> Bool
	Return a Boolean represented by string.
	"""
	if string.lower() == 'true':
		return True
	elif string.lower() == 'false' or string.lower() == '':
		return False
	else:
		raise ValueError("Cannot convert {} to a boolean".format(string))

# Data writing functions:
        
def write_data(data):
	rs = ("user,excess_locations,mean_distance,median_distance\n")
	for tup in data:
		user = tup[0]
		excess = round(tup[1], 2)
		mean = round(tup[2], 2)
		median = round(tup[3], 2)
		rs = rs + "{},{},{},{}\n".format(user,
		excess, mean, median)
	fd = open(output_compare_file, "w")
	fd.write(rs)

if __name__ == "__main__":
	location_data = compare_locations(locations_gt, output_locations_file)
	episode_data = compare_episodes(activities_gt, output_episodes_file)
	print(episode_data)
	write_data(location_data)
