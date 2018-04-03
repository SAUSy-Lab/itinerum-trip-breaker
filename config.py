#
# Configuration parameters 
#

# coordinates.csv file as from Itinerum
input_coordinates_file = './inputs/coordinates.csv'

# 
output_activities_file = "./outputs/activities.csv"
output_locations_file = "./outputs/locations.csv"
output_points_file = './outputs/points_diagnostic.csv'
time_file = "./outputs/location_times.csv"

# How much time must be spent in one spot for it to be detected as a potential  
# activity location? In seconds.
minimum_activity_time = 10*60

# Spatial kernel bandwidth in meters (standard deviation of gaussian kernel)
kernel_bandwidth = 25

# minimum distance between separate clusters 
# (parameter for activity location detection)
cluster_distance = 50

# what is the limit of stated h_accuracy which will be acceptable?
# (standard deviation in meters of a normal distribution?)
min_accuracy = 100

# interpolation distance parameter (meters). maximum length of segment to 
# remain uninterpolated for linear spatial interpolations
interpolation_distance = 30
