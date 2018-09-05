from points import Location
from config import *
from statistics import median
from datetime import timedelta, datetime
import editdistance


def read_headers(fname):
	""" (str) -> dict
	Return a dictionary mapping header names to column indices.
	Removes the need to hard code column numbers when reading files.
	
	See read_episodes for intended use.
	"""
	fd = open(fname)
	d = {}
	header = fd.readline()
	titles = header.split(",")
	for i in range(len(titles)):
		d[titles[i].strip()] = i
	fd.close()
	return d


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
			percent_excess = num_locations / len(true_locations[user])
			results.append((user, percent_excess, mean_min_dis, med))
	return results


def compare_user_locations(distances, num_to_match):
	""" (list, int) -> list
	Return a list of length num_to_match of the shortest distances
	between pairs of points represented by the list distances.

	This process, including parts of compare_locations where it is called,
	can be rewritenn to in O(n^2 log n) instead of O(n^3)
	using a recursive divide and conquer algorithm, standard for the closest pair problem. 
	"""
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
			# Remove any entries in which the computed or true location match
			if entry[1] == compd or entry[2] == truth:
				dist_list.remove(entry)
	return min_distances


def distance_matrix(h, matrix, truths, compds):
	""" (dict, dict, dict, dict) -> None
	Update matrix for user with a 2d matrix mapping true
	locations and computed locations to their respective distances.
	"""
	lid = h["location_id"]
	for location in truths:
		matrix[location[lid]] = {}
		for guess in compds:
			loc1 = Location( guess[h["lon"]], guess[h["lat"]] )
			loc2 = Location( location[h["lon"]], location[h["lat"]] )
			matrix[location[lid]][guess[lid]] = loc1.distance(loc2)


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
	""" (str, str) -> list
	Return a list of episode comparison metrics by user found in the 
	files truth and compd.

	A complete detailing of the metrics included is found
	in compare_user_episodes.
	"""
	true_users_to_eps = read_episodes(truth)
	computed_users_to_eps = read_episodes(compd)
	# get the headers, both files should match
	h = read_headers(truth)
	metrics = []
	for user in true_users_to_eps:
		if user in computed_users_to_eps:
			user_metrics = (user,) + compare_user_episodes(true_users_to_eps[user],
				computed_users_to_eps[user], h)
			metrics.append(user_metrics)
		else:
			print("{} not in computed episodes".format(user))
	return metrics


def read_episodes(file_name):
	""" (str) -> dict
	Return a dictionary mapping each user
	in file_name to a list of the episodes for that user.

	Episodes here are represented by a line from the episodes file,
	"""
	h = read_headers(file_name)
	fd = open(file_name)
	users_to_eps = {}
	fd.readline()
	for l in fd:
		line = l.strip().split(",")
		# use read_headers to avoid hard coding the column number
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


def compare_user_episodes(true, computed, h):
	""" (list, list, dict) -> tuple
	Return a tuple of episode comparison metrics for
	a particular user.

	The list values holds the intermediate values for
	directly computing the metrics.
	"""
	# TODO replace hardcoded unix time indices with header mapping
	end_time = min(true[-1][-1], computed[-1][-1])  # TODO
	start_time = max(true[0][-1], computed[0][-1])  # TODO
	total_time = float(end_time) - float(start_time)
	i, j = 0, 0
	# First we bring indices to overlapping episodes
	while true[i+1][-1] < start_time:  # TODO
		i = i + 1
	while computed[j+1][-1] < start_time:  # TODO
		j = j + 1

	# Values to compute:
	# 0 correct_unknown_time
	# 1 correct_travel_time
	# 2 correct_loc_time
	# 3 true_unknown_time
	# 4 true_travel_time
	# 5 true_loc_time
	# 6 incorr_unknown_time
	# 7 incorr_travel_time
	# 8 incorr_loc_time
	values = [0, 0, 0, 0, 0, 0, 0, 0, 0]
	# Iterate over episodes incrementally so that they always overlap
	# TODO minor bug when an episode has duration 0
	comp_str = ""
	true_str = ""
	while true[i][-1] < end_time and computed[j][-1] < end_time:
		duration = compare_single_episode((true[i], true[i+1]),
			(computed[j], computed[j+1]), h, values)
		if true[i+1][-1] <= computed[j+1][-1]:
			true_str = true_str + update_ep_str(true[i])
			i = i + 1
		else:
			comp_str = comp_str + update_ep_str(computed[j])
			j = j + 1
	# percent correct or incorrect unknown, travelling, or at_location time
	p_corr_ut = values[0] / values[3]  # percent correctly identified unknown time
	p_corr_trav = values[1] / values[4]  # percent correctly identified travelling time
	p_corr_loc = values[2] / values[5]  # percent correctly identified time at a location
	p_inc_ut = values[6] / (total_time - values[0])  # percent incorrectly identified unknown time
	p_inc_trav = values[7] / (total_time - values[1])  # percent incorrectly identified travelling time
	p_inc_loc = values[8] / (total_time - values[2])  # percent incorrectly identified time spent at location
	return (p_corr_ut, p_corr_trav, p_corr_loc, p_inc_ut, p_inc_trav,
		p_inc_loc, true_str, comp_str)


def compare_single_episode(true_pair, computed_pair, h, values):
	"""
	Update the appropriate episode comparison metrics in values,
	and return the amount of time the episodes true_pair and computed_pair
	overlap.
	"""
	overlapping_time = float(min(true_pair[1][-1],
		computed_pair[1][-1])) - float(max(true_pair[0][-1], computed_pair[0][-1]))
	true_unknown = False if true_pair[0][h["unknown"]] == "" else True
	true_at_loc = False if (true_pair[0][h["location_id"]] == "" or
		true_pair[0][h["location_id"]] == "None") else True
	comp_unknown = False if computed_pair[0][h["unknown"]] == "" else True
	comp_at_loc = False if (computed_pair[0][h["location_id"]] == "" or
		computed_pair[0][h["location_id"]] == "None")else True

	if true_unknown and comp_unknown:
		values[0] += overlapping_time
	if not (true_at_loc or comp_at_loc):
		values[1] += overlapping_time
	if true_at_loc and comp_at_loc:
		values[2] += overlapping_time
	if true_unknown:
		values[3] += overlapping_time
	if not true_at_loc:
		values[4] += overlapping_time
	else:
		values[5] += overlapping_time
	if not true_unknown and comp_unknown:
		values[6] += overlapping_time
	if true_at_loc and not comp_at_loc:
		values[7] += overlapping_time
	if not true_at_loc and comp_at_loc:
		values[8] += overlapping_time

	return overlapping_time


def update_ep_str(episode):
	"""
	Return a letter representing whether
	this episode is classified as unknown time,
	travelling time, or activity time (time spent at a location)
	"""
	if episode[4] == "True":
		return "U"
	elif episode[2] != "":
		return "A"
	else:  # travel
		return "T"

# Data writing functions:


def write_data(data):
	""" (tuple) -> NoneType
	Write out the contents of data to the .csv file
	specified in config.py
	"""
	rs = ("user,percent_excess_locations,mean_distance,median_distance," +
		"percent_identified_unknown,percent_identified_travel," +
		"percent_identified_location,percent_misidentified_unknown," +
		"percent_misidentified_travel,percent_misidentified_location," +
		"true_eps,comp_eps,edit_distance,record_length,computed_length\n")
	for tup in data:
		user = tup[0]
		excess = round(tup[1], 2)
		mean = round(tup[2], 2)
		median = round(tup[3], 2)
		p_i_ut = round(tup[4], 2)
		p_i_tt = round(tup[5], 2)
		p_i_lt = round(tup[6], 2)
		p_m_ut = round(tup[7], 2)
		p_m_tt = round(tup[8], 2)
		p_m_lt = round(tup[9], 2)
		comp_string = tup[11]
		true_string = tup[10]
		edit_distance = editdistance.eval(comp_string, true_string)
		record_len = len(true_string)
		computed_len = len(comp_string)
		rs = rs + "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(user,
			excess, mean, median, p_i_ut, p_i_tt, p_i_lt, p_m_ut, p_m_tt,
			p_m_lt, true_string, comp_string, edit_distance, record_len, computed_len)
	fd = open(output_dir+'/compare.csv', "w")
	fd.write(rs)


def merge_lists(t1, t2):
	""" (tuple, typle) -> List
	Merge two lists, omitting the user_id from the second
	in order to facilitate writing to a file.
	"""
	new_data = []
	for i in range(len(t1)):
		new_data.append(t1[i] + t2[i][1:])
	return new_data

if __name__ == "__main__":
	# Collect location comparison metrics by user
	location_data = sorted(compare_locations(locations_gt,
		output_dir+'/locations.csv'))
	# Collect episode comparison metrics by user
	episode_data = sorted(compare_episodes(activities_gt,
		output_dir+'/episodes.csv'))
	# Merge collected data and write out.
	data = merge_lists(location_data, episode_data)
	write_data(data)
