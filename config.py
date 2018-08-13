#
# Configuration parameters
#

# Flag for toggling debugging print statement on or off
db_out = True

input_dir = './inputs'
output_dir = './outputs'

# input data directly from Itinerum
input_coordinates_file = input_dir + '/coordinates.csv'
input_survey_responses_file = input_dir + '/survey_responses.csv'

# output files generated by these scripts
output_episodes_file = output_dir + '/episodes.csv'
output_locations_file = output_dir + '/locations.csv'
output_points_file = output_dir + '/classified_points.csv'
output_days_file = output_dir + '/days.csv'
output_compare_file = output_dir + '/compare.csv'

# manual ground truth data for comparison with compare.py
locations_gt = output_dir + '/locations_ground_truth.csv'
activities_gt = output_dir + '/episodes_ground_truth.csv'

# timezone to use (until efficient timezone lookups can be implemented)
local_timezone = 'America/Toronto'

# How much time must be spent in one spot for it to be detected as 
# an activty episode?
minimum_activity_time = 10*60 # seconds

# Spatial kernel bandwidth in meters (standard deviation of gaussian kernel)
# Used for location detection.
kernel_bandwidth = 25

# Limit of h_accuracy beyond which points get discarded
# (is the unit  meters of in standard deviation of a normal distribution?)
min_accuracy = 100

# Minimum distance between separate locations
location_distance = 150  # meters

# interpolation distance parameter (meters). maximum length of segment to
# remain uninterpolated for linear spatial interpolations.
# For reasonable results, this must be < cluster_distance
interpolation_distance = 30
assert location_distance > interpolation_distance

# weight coefficient:
# Affects the input to the weight distribution function,
# which varies from 0.5 to 1 as the ratio of distance over time difference
# increases. Increasing this value over 1 scales up this ratio
# while decreasing under 1 scales increases it.
weight_coef = 1
assert (weight_coef > 0)

# Flag for toggling different denominators in calculating episode detection
# Metrics. If on, the denominator consists of the total time for that user,
# If off, it consists of only the portions of episodes that could have affected
# that particular metric, in particular:
# correctly identified unknown time and activity time (not location dependent)
# and misidentified unknown and activity time.
percent_total = True

# Number of worker processes on which to run main.py
# and a flag toggling whether or not to use multiprocessing for main.py
from os import cpu_count
num_pro = cpu_count()
multi_process = False
