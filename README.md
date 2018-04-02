## Introduction
This is an application for parsing [Itinerum](https://github.com/TRIP-Lab/itinerum-android) travel survey data into a more standard travel survey format, showing trips, activity locations, modes of travel and times of arrival/departure, etc. 

## Algorithm
For each user:

### Data cleaning
The basic idea in the cleaning phase to remove points not based on decent GPS data. Many points may be derived from cell-towers, wifi networks, etc. The phone is a black box in this regard. Such points often repeat a location precisely, which is unlikely with a GPS reading or they may suddenly appear very far away from their temporal neighbors. 

1. Remove points with high known error (h_accuracy > x meters) 
2. Remove points at same location as temporal neighbors
3. Any major jump away and back again, especially if the h_error doesn't justify the distance may indicate a non-GPS signal or a bad error estimate. Away-and-back-again points are identified by the minimum distance from neighboring points and the angle formed between the three. 
4. Repeat step 2
5. More to come...

### Data segmentation
Data is segmented into lists of points where we feel confident that we haven't lost track of the user, e.g. through a powered-off phone. We refer to these as *known segments* and let the surrounding time be considered *unknown*. For the moment, this consists only of:

1. Consider as unknown any time where the user moves more than 1 kilometer and two hours without reporting a location. 

This is obviously not sufficient, and there is much more to come here. For example, subways will need to be scrutinized. 

### Location detection
This phase largely follows the [method described by Thierry Chaix and Kestens](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3637118/). We essentially do a time-weighted KDE on the user's points, spatio-temporally linear-interpolated where necessary. The kernel density estimation will be done in R (for lack of a decent Python package) and then brought back into Python. 

1. Calculate a KDE based on time weighted GPS points, interpolated as appropriate. GPS points are also used as points where the PDF is estimated. This saves the cost of estimating a grid of points over sparse GPS data. 
2. Estimate an "activity threshold" for that surface. Id est, how high would a peak in the KDE be if someone spent X seconds at a given location, given the known parameters of 
    - total time under the surface (sum of time weights)
    - average GPS error for the user
    - kernel bandwidth
3. Points with a PDF estimate above the threshold are clustered into contiguous groups.
4. The maximum PDF estimate from each cluster (the peak) will be taken as a potential activity location. (It's possible that it may make sense to use a polygonized version of the cluster as a definition of the activity location rather than a point from the peak. This is not implemented yet.)

### Activity/Trip sequencing
The ultimate goal of this program is to create a travel-diary-like sequence of activities; this step is where that happens. Conceptually, time may be spent in only one of two (three) ways in this model: 1) travelling, 2) at an 'activity' (or 3) it may be classified as unknown). Framing the problem thus leads to a [hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) approach to discrete activity sequence classification. Each observed GPS point can be assigned an emission probability associated with each location, based on distance between the two. Transition probabilities give strong preference to state continuity and disallow teleporting (activities without intermediating travel).

The Viterbi algorithm provides the most likely underlying classification of points, and these are used to map time to travel and activities. 


## Dependencies
* Python 3
    - Rpy2 
    - pyproj
    - scipy
* R 3.3+
    - ks package
