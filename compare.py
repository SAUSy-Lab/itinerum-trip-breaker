from point import Point
from config import *
from misc_funcs import distance
from statistics import median
from datetime import timedelta, datetime
import misc_functions


# TODO this should be refactored
def compare_locations(truth, compd):
	""" (str, str) -> [(str, int, float, float)]
	Compare ground truth locations to the computed locations,
	and output a list of quality metrics by user.
	"""
	true_locations = {}
	computed_locations = {}
	get_locations(truth, true_locations)
	get_locations(compd, computed_locations)

	users_to_matrix = {}
	results = []
	for user in true_locations:
		if user not in computed_locations:
			print(user, " not in computed locations")
		else:
			num_locations = len(computed_locations[user]) - len(true_locations[user])
			users_to_matrix[user] = {}
			distance_matrix(
                            user, users_to_matrix[user], true_locations[user], computed_locations[user])
				# users_to_matrix[user][true_location][computed_location] = distance
				# dictionary of dictionaries of dictionaries of distances
			min_distances = []
			dist_list = []
			for loc in users_to_matrix[user]:
				for guess in users_to_matrix[user][loc]:
					dist_list.append((users_to_matrix[user][loc][guess], guess, loc))
			while remaining_locations(dist_list, num_locations):
				dist_list.sort()
				best = dist_list.pop(0)
				min_distances.append(best[0])
				guess = best[1]
				loc = best[2]
				for entry in dist_list:
					if entry[1] == guess or entry[2] == loc:
						dist_list.remove(entry)
			mean_min_dis = sum(min_distances) / len(min_distances)
			med = median(min_distances)
			results.append((user, num_locations, mean_min_dis, med))
	return results

def remaining_locations(dist_list, excess):
	""" ([(float, int, int)]) -> Bool
	Return true iff there are more remaining unique locations in dist_list than
	we expect, given our original excess locations.
	"""
	locs = []
	guess = []
	for entry in dist_list:
		if entry[1] not in guess:
			guess.append(entry[1])
		if entry[2] not in locs:
			locs.append(entry[2])
	if excess > 0:
		return len(guess) > excess
	elif excess < 0:
		return len(locs) < excess
	else: #excess ==0
		return len(locs) == len(guess) == 0
	
def distance_matrix(user, matrix, truths, compds):
	""" (str, dict, dict, dict) -> None
	Update matrix for user with a 2d matrix mapping true
	locations and computed locations to their respective distances.
	"""
	for locus in truths:
		matrix[locus[1]] = {}
		for guess in compds:
			# TODO data fields shouldn't be hardcoded
			p1 = Point("", guess[2], guess[3], 0)
			p2 = Point("", locus[2], locus[3], 0)
			# Felipevh forgets if these need to be projected first
			matrix[locus[1]][guess[1]] = distance(p1, p2) 

def get_locations(location_file, utl):
	""" (str, dict) -> None
	Populate utl with the locations described in location_file.
	"""
	fd = open(location_file)
	line = fd.readline().split(',')
	for line in fd:
		loc = line.split(',')
		if loc[0] not in utl:
			utl[loc[0]] = [loc]
		else:
			utl[loc[0]].append(loc)
                    
def compare_episodes(truth, guess):
	""" (str, str) -> []
	"""
	truth_dict = read_file(truth)
	guess_dict = read_file(guess)
	# True and computed unknown times
	tut = find_unknown_time(truth_dict[5], truth_dict[4]) # TODO don't hardcode
	cut = find_unknown_time(gues_dict[5], guess_dict[4])

def read_file(fname):
	""" (str) -> {int : [str]}
	Return a dictionary mapping column numbers
	to lists of entries for that column in fname.
	Drops the header line and strips whitespace.
	"""
	pass

def find_unknown_time(start_times, uflags):
	""" ([str], [str]) -> [(timedelta, bool)]
	Return a list of timedeltas, and true iff
	that time is classified as unknown
	"""
	assert(len(start_time) == len(uflags))
        result = []
	for i in range(len(start_times)):
		if start_times[i].endswith("M"):
			ts = parse_gt_ts(start_times[i])
		else:
			ts, _ = parse_ts(start_times[i] + "-00:00")
		result.append((ts, literal_eval(uflags[i])))
	return result

def parse_gt_ts(t):
	""" (str) -> DateTime
	'mm/dd/yyyy HH:mm XX'
	'yyyy-mm-ddTHH:mm:ss-00:00'

	Parse archaic timestamps into a human readable format.
	"""
	real_timestamp = t[6:10]+"-"+t[:2]+"-"t[3:5]+"T"+t[11:13]+":"+t[14:16]+"00-00:00"
	return misc_functions(real_timestamp)[0]

def literal_eval(string):
	""" (str) -> Bool
	Return a Boolean represented by string.
	"""
	if string.lower() == 'true':
		return True:
	elif string.lower() == 'false' or string.lower() == '':
		return False
	else:
		raise ValueError("Cannot convert {} to a boolean".format(string))

if __name__ == "__main__":
	print(compare_locations(locations_gt, output_locations_file))
	#print(compare_episodes(config.activities_gt, config.output_episodes_file))
