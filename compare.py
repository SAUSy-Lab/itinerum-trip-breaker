from point import Point
from config import *
from misc_funcs import distance
from statistics import median
from datetime import timedelta, datetime
from misc_funcs import parse_ts, read_headers


def compare_user_locations(distances, num_to_match):
	min_distances = []
	dist_list = []
	for loc in distances:
		for guess in distances[loc]:
			dist_list.append((distances[loc][guess], guess, loc))
	for _ in range(num_to_match):
		dist_list.sort()
		best = dist_list.pop(0)  # First in the list
		min_distances.append(best[0])  # the distance entry
		guess = best[1]
		loc = best[2]
		for entry in dist_list:
			# if the guess or true location match
			if entry[1] == guess or entry[2] == loc:
				dist_list.remove(entry)
	return min_distances


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


def distance_matrix(h, matrix, truths, compds):
	""" (dict, dict, dict, dict) -> None
	Update matrix for user with a 2d matrix mapping true
	locations and computed locations to their respective distances.
	"""
	for locus in truths:
		matrix[locus[1]] = {}  # what is this 1 value? TODO
		for guess in compds:
			p1 = Point("", guess[h["lon"]], guess[h["lat"]], 0)
			p2 = Point("", locus[h["lon"]], locus[h["lat"]], 0)
			# Felipevh forgets if these need to be projected first
			matrix[locus[1]][guess[1]] = distance(p1, p2)


def get_locations(location_file, utl):
	""" (str, dict) -> None
	Populate utl with the locations described in location_file.
	"""
	fd = open(location_file)
	h = read_headers(location_file)
	line = fd.readline().split(',')
	for line in fd:
		loc = line.split(',')
		if loc[h["user_id"]] not in utl:
			utl[loc[h["user_id"]]] = [loc]
		else:
			utl[loc[h["user_id"]]].append(loc)


def compare_episode_time(truth, guess):
	""" (str, str) -> [str, float, float]
	Return a list of users and their episode time quality metrics.
	"""
	result = []
	truth_dict = read_file(truth)
	guess_dict = read_file(guess)
	# True and computed unknown times
	ht = read_headers(truth)
	hg = read_headers(guess)
	tut = find_unknown_time(truth_dict[ht["start_time"]],
				truth_dict[ht["unknown"]], truth_dict[ht["user_id"]])
	cut = find_unknown_time(guess_dict[hg["start_time"]],
				guess_dict[hg["unknown"]], guess_dict[hg["user_id"]])
	for user in tut.keys():
		if user not in cut.keys():
			print("user {} not in computed data".format(user))
		else:
			result.append((user, compare_user_eps(tut[user], cut[user])))
	return result


def compare_user_eps(truth, computed):
	""" ([(Datetime, Bool)], [(Datetime, Bool)]) -> (float, float)
	Return the episode time quality metrics for this user.
	"""
	start_time = max(truth[0][0], computed[0][0])
	end_time = min(truth[-1][0], computed[-1][0])
	# Time in minutes that the survey lasted, trimmed
	##total_time = ((end_time - start_time).days * 24 * 60 +
	##(end_time - start_time).seconds / 60)
	total_time = 0
	i, j = 0, 0
	# total unknown time, correctly identified, misidentified
	tut, ciut, mut = 0, 0, 0
	# activity time
	tat, ciat, mat = 1, 0, 0
	# while both iterators haven't reached the end of the survey time
	while truth[i][0] < end_time and computed[j][0] < end_time:
		if overlaps(truth[i][0], truth[i+1][0], computed[j][0], computed[j+1][0]):
			# identified time lies between latest begining and earliest end
			ee = min(truth[i+1][0], computed[j+1][0])
			lb = max(truth[i][0], computed[j][0])
			it = (ee - lb)
			if truth[i][1] and computed[j][1]:  # correctly identified unknown time
				ciut += it.days * 24 * 60 + it.seconds / 60
			if not truth[i][1] and computed[j][1]:  # misidentified unknown time
				mut += it.days * 24 * 60 + it.seconds / 60
			if not (truth[i][1] or computed[j][1]):  # correctly identified activity time
				ciat += it.days * 24 * 60 + it.seconds / 60
			if truth[i][1] and not computed[j][1]:  # misidentified unknown time
				mat += it.days * 24 * 60 + it.seconds / 60
		# whichever episode ends first needs to be incremented
		if truth[i+1][0] < computed[j+1][0]:  # truth ep ends first
			t = ((truth[i+1][0] - truth[i][0]).days * 24 * 60 +
			     (truth[i+1][0] - truth[i][0]).seconds / 60)
			if truth[i][1]:
				tut += t
				total_time += t
			else:
				tat += t
				total_time += t
			i += 1
		else:  # computed ep ends first
			j += 1
	pciut = ciut / tut
	pmiut = mut / (total_time - tut)
	print("ciat: {}, tat: {}, mat: {}, tot: {}".format(ciat, tat, mat, total_time))
	pciat = ciat / tat
	pmat = mat / (total_time - tat)
	return (pciut * 100, pmiut * 100, pciat * 100, pmat * 100)


def overlaps(t1, t2, c1, c2):
	""" (Datetime, Datetime, Datetime, Datetime) -> Bool
	Return true iff the time gaps t2 - t1 and c2 - c1 overlap.
	"""
	return (t1 < c1 and t2 > c1) or (c1 < t1 and c2 > t1)


def read_file(fname):
	""" (str) -> {int : [str]}
	Return a dictionary mapping column numbers
	to lists of entries for that column in fname.
	Drops the header line and strips whitespace.
	"""
	fd = open(fname, "r")
	fd.readline()  # header
	d = {}
	for line in fd:
		cleaned = line.split(',')
		for i in range(len(cleaned)):
			if i in d.keys():
				d[i].append(cleaned[i].strip())
			else:
				d[i] = [cleaned[i].strip()]
	return d


def find_unknown_time(start_times, uflags, users):
	""" ([str], [str]) -> [(timedelta, bool)]
	Return a list of timedeltas, and true iff
	that time is classified as unknown
	"""
	assert(len(start_times) == len(uflags))
	result = {}
	for i in range(len(start_times)):
		lst = []
		user = users[i]
		if start_times[i].strip().endswith("M"):
			ts = parse_gt_ts(start_times[i])
			lst.append((ts, literal_eval(uflags[i])))
		elif start_times[i] == "":
			pass
		else:
			ts, _ = parse_ts(start_times[i] + "-00:00")
			lst.append((ts, literal_eval(uflags[i])))
		if user not in result.keys():
			result[user] = lst
		else:
			result[user].extend(lst)
	return result


# TODO remove this once timestamps are standardized to epoch time
def parse_gt_ts(t):
	""" (str) -> DateTime
	'mm/dd/yyyy HH:mm XX'
	'yyyy-mm-ddTHH:mm:ss-00:00'

	Parse archaic timestamps into a human readable format.
	"""
	year = t[6:8]
	month = t[:2]
	day = t[3:5]
	minute = t[12:14]
	hour = t[9:11]
	second = "01"
	real_timestamp = "20{}-{}-{}T{}:{}:{}-00:00".format(year, month,
							day, hour, minute, second)
	return parse_ts(real_timestamp)[0]


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


def write_data(locs, eps):
	user_to_loc = {user: [excess, mean, med]
		for (user, excess, mean, med) in locs}
	user_to_eps = {user: [ciut, mut, ciat, mat] for (user, (ciut, mut, ciat, mat)) in eps}
	rs = ("user,excess_locations,mean_distance,median distance" +
		",identified_unknowntime,misidentified_unknowntime" +
		",identified_activitytime,misidentified_activitytime\n")
	for user in user_to_loc.keys():
		excess = user_to_loc[user][0]
		mean = round(user_to_loc[user][1], 2)
		median = round(user_to_loc[user][2], 2)
		ciut = round(user_to_eps[user][0], 2)
		mut = round(user_to_eps[user][1], 2)
		ciat = round(user_to_eps[user][2], 2)
		mat = round(user_to_eps[user][3], 2)
		rs = rs + "{},{},{},{},{},{},{},{}\n".format(user,
		excess, mean, median, ciut, mut, ciat, mat)
	fd = open(output_compare_file, "w")
	fd.write(rs)

if __name__ == "__main__":
	loc = compare_locations(locations_gt, output_locations_file)
	eps = compare_episode_time(activities_gt, output_episodes_file)
	write_data(loc, eps)
