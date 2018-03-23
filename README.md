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

1. Calculate a KDE based on time weighted GPS points, interpolated as appropriate. GPS points are also used as points where the PDF is estimated. This saves the cost of estimating a grid of points over sparse GPS data. 
2. Estimate an "activity threshold" for that surface. Id est, how high would a peak in the KDE be if someone spent X seconds at a given location, given the known parameters of 
    - total time under the surface (sum of time weights)
    - average GPS error for the user
    - kernel bandwidth
3. Points with a PDF estimate above the threshold are clustered into contiguous groups.
4. The maximum PDF estimate from each group (the peak) will be taken as a possible activity location. (It's possible that it may make sense to use a polygonized version of the cluster as a definition of the activity location rather than a point from the peak. This is not implemented yet.)

### Activity/Trip sequencing
This phase is a bit hazier as it's still two steps away from implementation, but the idea is as follows:
* Treat time sufficiently near a potential activity location as an activity
* Treat time between activity locations as a trip
* Look for activity/travel times that are implausibly long/short and reallocate time if necessary.
* Potential activity locations that never get visited or get visited for many brief periods only will be ignored and time reallocated to trips, etc.


## Dependencies
* Python 3
    - Rpy2 (python module)
    - pyproj (python module)
    - scipy (python module)
* R 3.3+
    - ks (R package)
