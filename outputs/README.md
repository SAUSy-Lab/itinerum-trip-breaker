This data comes from one week of a SAUSy researcher using the Itinerum app and is provided for testing purposes only.

*Please be very careful not to ever share any confidential data from survey participants here or elsewhere.*

## Output Files
Output files generated will include:

### locations.csv
This file gives activity locations for each user. Activity locations are points. 

It includes the following feilds: `user_id,location_id,lon,lat,description,used`.
* `used` indicates whether the activity location was actually visited. It's possible to detect an activity location but not identify any activites there if there are repeated short duration visits there, e.g a bus stop. 


### episodes.csv
This file gives the sequence of travel, activities, and unknown times as the user proceeds through time. Each entry is the *start* of an episode. Each user appears to us out of unknown time and then begins some series of activities before ending again by disapearing into unknown time. Unknown is a binary classification showing whether we think we know what the user is doing. A user in unknown time may have some data points but not enough to infer anything. 

Fields include `user_id,sequence,location_id,mode,unknown,start_time`
* `sequence` increments per user.
* `location_id` refers to `locations.csv`. It's not decided yet whether location_id is incremented per user or globally.
* `mode` is not currently used except in our ground truth data. The idea is that this will hold a delimited string indicating a series of modes of travel, e.g. walk/transit/walk.
* `start_time` the format for this field is not yet set in stone. 


### days.csv
This file summarizes `episodes.csv` per day in the dataset. We treat days not as an exact clock day, but extend the day a few hours past midnight. We do this on the assumption that some people stay up past midnight, but very few people will be up and moving about at 3 or 4 AM. 

Fields include: `user_id,date,DoW,total_minutes,trip_count,travel_time,unknown_time,home_time,work_time,school_time,home_count,work_count,school_count`
* `DoW` is an integer indicating the day of the week, with Sunday = 0. 
* The count fields are counts of activities at the known location types. 
* The time feilds are given in minutes and should sum to 1440, the number of minutes in a day, except for the first and last day in the dataset. 


### classified_points.csv
This file is meant for debugging output and eventually for returning trip geometry. It includes the fields `user_id,lon,lat,weight,removed,interpolated,state`. 
* `removed` means the point was dropped during the cleaning stage
* `interpolated` means that the point was generated synthetically and was not in the original input data
* `state` is the hidden Markov classification assigned to the point. This is an integer feild referring to activity location id's for a user except for `0` which is travel. 
