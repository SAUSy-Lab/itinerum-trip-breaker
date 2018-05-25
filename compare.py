from point import Point
import config
from misc_funcs import distance
# what does activities_gt get compared to? TODO

def compare_locations(truth, compd):
	true_locations = {}
	computed_locations = {}
	get_locations(config.locations_gt, true_locations)
	get_locations(compd, computed_locations)

	users_to_matrix = {}
	for user in true_locations:
		if user not in computed_locations:
			print(user, " not in computed locations")
		else:
			size = max(len(true_locations[user]), len(computed_locations[user]))
			users_to_matrix[user] = {}
			distance_matrix(
                            user, users_to_matrix[user], true_locations[user], computed_locations[user])
				# users_to_matrix[user][true_location][computed_location] = distance
				# dictionary of dictionaries of dictionaries of distances

	print(users_to_matrix)
def distance_matrix(user, matrix, truths, compds):
	for loci in truths:
		matrix[loci[1]] = {}
		for guess in compds:
			# TODO data fields shouldn't be hardcoded
			p1 = Point(ts, guess[2], guess[3], 0)
			p2 = Point("", loci[2], loci[3], 0)
			# Felipe forgets if these need to be projected first
			matrix[loci[1]][guess[1]] = distance(p1, p2) 

def get_locations(location_file, utl):
	fd = open(location_file)
	line = fd.readline().split('\n')
	for line in fd:
		loc = line.split('\n')
		if loc[0] not in utl:
			utl[line[0]] = [loc]
		else:
			utl[line[0]].append(loc)
                    

if __name__ == "__main__":
	compare_locations(config.locations_gt, config.output_locations_file)    
