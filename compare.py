from point import Point
import config
from misc_funcs import distance
# what does activities_gt get compared to? TODO

def compare_locations(truth, compd):
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
				best = dist_list.pop(0) #removes all the ground truth truth entries
				min_distances.append(best[0])
				guess = best[1]
				loc = best[2]
				for entry in dist_list:
					if entry[1] == guess or entry[2] == loc:
						dist_list.remove(entry)
			min_avg_dis = sum(min_distances) / len(min_distances)
			median = 0 #TODO
			results.append((user, num_locations, min_avg_dis, median))
	return results

def remaining_locations(dist_list, excess):
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
	for locus in truths:
		matrix[locus[1]] = {}
		for guess in compds:
			# TODO data fields shouldn't be hardcoded
			p1 = Point("", guess[2], guess[3], 0)
			p2 = Point("", locus[2], locus[3], 0)
			# Felipe forgets if these need to be projected first
			matrix[locus[1]][guess[1]] = distance(p1, p2) 

def get_locations(location_file, utl):
	fd = open(location_file)
	line = fd.readline().split(',')
	for line in fd:
		loc = line.split(',')
		if loc[0] not in utl:
			utl[loc[0]] = [loc]
		else:
			utl[loc[0]].append(loc)
                    

if __name__ == "__main__":
	print(compare_locations(config.locations_gt, config.output_locations_file))    
