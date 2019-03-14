This data comes from one week of a SAUSy researcher using the Itinerum app and is provided for testing purposes only.

*Please be very careful not to ever share any confidential data from survey participants here or elsewhere.*
*As a precaution, git will ignore all .csv files in the outputs folder, other than those already present for testing*
## Output Files
Output files generated will include:

### locations.csv
This file gives activity locations for each user. Activity locations are points. 

It includes the following feilds: `user_id,location_id,lon,lat,description,used`.
* `used` indicates whether the activity location was actually visited. It's possible to detect an activity location but not identify any activites there if there are repeated short duration visits there, e.g a bus stop. 
* `description` may contain one of `('home','work','school')` if named locations were provided per user. Otherwise this will be blank.


### episodes.csv
This file gives the sequence of travel, activities, and unknown times (collectively 'episodes') as the user proceeds through time. Each entry is the *start* of an episode. Each user appears to us by beginning a travel or activity episode followed by one or more additional episodes before disapearing into a terminal 'unknown' episode. Unknown is a binary indicating whether or not we think we know what the user is doing. A user in unknown time may have some data points but not enough to infer anything. Or they have finished the survey and carried on with their life. In this case, the timestamp just indicates the end of the previous episode. 

Fields include `user_id,sequence,location_id,mode,unknown,local_start_time,unix_start_time,duration`
* `sequence` simply increments per user, helping to order records
* `location_id` refers to `locations.csv`. Id is incremented per user and is not globally unique.
* `mode` is not currently used except in our ground truth data. The idea is that this will hold a delimited string indicating a series of modes of travel, e.g. walk/transit/walk.
* `local_start_time` local timestamp with timezone (e.g. `2016-09-02 16:44:55-04:00`) indicating the *start* time of an episode.
* `unix_start_time` same time as the above, but formatted as a unix timestamp (e.g. `1472849095`)
* `duration` duration of the episode in seconds. 


### days.csv
This file summarizes `episodes.csv` per day in the dataset. We treat days not as an exact clock day, but extend the day a few hours past midnight. We do this on the assumption that some people stay up past midnight, but very few people will be up and moving about at 3 or 4 AM. 

Fields include: `user_id,date,DoW,total_minutes,trip_count,travel_time,unknown_time,home_time,work_time,school_time,home_count,work_count,school_count`
* `DoW` is an integer indicating the day of the week, with Sunday = 0, Monday = 1
* The count fields are counts of activities at the known location types. 
* The time fields are given in minutes and should sum to 1440, the number of minutes in a day, except for the first and last day in the dataset. This may vary during clock changes as these are based on local clock time. 


### classified_points.csv
This file is meant for debugging output and eventually for returning trip geometry. It includes the fields `user_id,lon,lat,weight,removed,interpolated,state`. 
* `removed` means the point was dropped during the cleaning stage
* `interpolated` means that the point was generated synthetically and was not in the original input data
* `state` is the hidden Markov classification assigned to the point. This is an integer feild referring to activity location id's for a user except for `0` which is travel. 
