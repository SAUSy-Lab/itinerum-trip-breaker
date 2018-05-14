#
# Configuration parameters 
#

# coordinates.csv file as from Itinerum
input_coordinates_file = './inputs/coordinates.csv'
input_survey_responses_file = './inputs/survey_responses.csv'

# sequential episodes file 
output_episodes_file = "./outputs/episodes.csv"
# locations file
output_locations_file = "./outputs/locations.csv"
# diagnostic output
output_points_file = './outputs/classified_points.csv'
# day by day summary file
output_days_file = "./outputs/days.csv"

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
interpolation_distance = 50
