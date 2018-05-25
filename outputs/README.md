This data comes from one week of a SAUSy researcher using the Itinerum app and is provided for testing purposes only.

*Please be very careful not to ever share any confidential data from survey participants here or elsewhere.*

## Output Files
Output files generated will include:

### locations.csv

### days.csv

### episode.csv

### classified_points.csv
This file is meant for debugging output and eventually for returning trip geometry. It includes the fields `user_id,lon,lat,weight,removed,interpolated,state`. 
* `removed` means the point was dropped during the cleaning stage
* `interpolated` means that the point was generated synthetically and was not in the original input data
* `state` is the hidden Markov classification assigned to the point. This is an integer feild referring to activity location id's for a user except for `0` which is travel. 
