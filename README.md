## Introduction
This is an application for parsing Itinerum travel survey data into a more standard travel survey format, showing trips, activity locations, modes of travel and times of arrival/departure, etc. 

## Algorithm
For each user:

### Data cleaning
The basic idea in the cleaning phase to remove points not based on decent GPS data. Many points may be derived from cell-towers, wifi networks, etc. The phone is a black box in this regard. Such points often repeat a location precisely, which is unlikely with a GPS reading or they may suddenly appear very far away from their temporal neighbors. 

1. Remove points with high known error (h_accuracy > x meters) 
2. Remove points at same location as temporal neighbors
3. Remove points based on angle and distance to temporal neighbors
4. Repeat step 2
5. More to come...

### Data segmentation
Data is segmented into lists of points where we feel confident that we haven't lost track of the user, e.g. through a powered-off phone. We refer to these as *known segments* and let the surrounding time be considered *unknown*. For the moment, this consists only of:

1. Consider as unknown any time where the user moves more than 1 kilometer and two hours without reporting a location. 

This is obviously not sufficient, and there is much more to come here. 

### Location detection
This phase largely follows the [method described by Thierry Chaix and Kestens](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3637118/). We essentially do a time-weighted KDE on the user's points, spatio-temporally linear-interpolated where necessary. The kernel density estimation will be done in R (for lack of a decent Python package) and then brought back into Python. 

The plan (not yet implemented) is to:
1. Scale the area under the probability surface to equal the time in the user's known segments. 
2. Set a threshold value based on cell-size to determine where sufficient time has been spent to suggest an activity.
3. Cells meeting the threshold will be clustered into contiguous groups.
4. The maximum of these groups (peaks) will be taken as a possible activity location. (It's possible that it may make sense to use a polygonized version of the cluster as a definition of the activity location rather than a point from the peak.)

### Activity/Trip sequencing
This phase is a bit hazier as it's still two steps away from implementation, but the idea is as follows:
* Treat time sufficiently near a potential activity location as an activity
* Treat time between activity locations as a trip
* Look for activity/travel times that are implausibly long/short and reallocate time if necessary.
* Potential activity locations that never get visited or get visited for many brief periods only will be ignored and time reallocated to trips, etc.


## Dependencies
Python 3

Rpy2 (python module)
pyproj (python module)
scipy (python module)

R 3.3+
ks (R package)
