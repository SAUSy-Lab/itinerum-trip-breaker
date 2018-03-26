#
# Configuration parameters 
#

# coordinates.csv file as from Itinerum
input_coordinates_file = 'inputs/coordinates.csv'

# 
FILENAME = "./outputs/activities.csv"
FILENAME_L = "./outputs/locations.csv"

# How much time must be spent in one spot for it to be detected as a potential  
# activity location? In seconds.
MIN_SECS_AT_LOC = 10*60

# Spatial kernel bandwidth in meters (standard deviation of gaussian kernel)
BANDWIDTH = 25

# minimum distance between separate clusters 
# (parameter for activity location detection)
CLUSTER_DISTANCE = 50
