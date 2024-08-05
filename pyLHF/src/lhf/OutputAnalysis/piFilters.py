#!/usr/bin/python3

from kneed import KneeLocator
import numpy as np
import bisect
import sys

#### Set of functions to filter sets of persistence intervals (PIs) of the form: <dim, birth, death>
#### The filtering will be performed based on the length (<death> - <birth>) of the PIs.  The filtering can be performed uniformly
#### across all dimensions of the PIs or it can be directed to perform the filtering uniquely to the PIs at each dimension.

def _filterPIs(data, filter=('elbow', 1.0)) :
    """
    this is a private function that we use to actually filter the data.  this function actually implements all of the filtering
    operations in one common code base.  in this case, the dimension field of the PIs is ignored, so to filter by dimension, this
    function must be called with only the data to be filtered (thus, you can filter all dimensions, one dimension, or a set of
    dimensions).   

    Parameters
    ----------
    data : numpy.ndarray
        numpy [*, 3] matrix; each row is a persistence interval written as <dimension, birth, death>.
    filter : tuple: (filterType:str , param:float)  
        type of filter to apply, see details in public function 'filterPIs'

    Returns
    -------
    filtered persistence intervals
    """

    def applyCutoff (data, cutoff) :
        if cutoff > data.shape[0] or cutoff < 0 :
            cutoff = data.shape[0]
        return data[cutoff:, :]

    # don't filter PI sets that have less than 1 element.
    if not len(data) > 1 :
        return data

    # process the filter argument string/list
    filterType = None
    filterArg1 = None
    filterArg2 = None
    if isinstance(filter, str) :
        filterType = filter
    else :
        filterType = filter[0]
        filterArg1 = filter[1]
        if len(filter) == 3 :
            filterArg = filter[2]

    # sort the persistence intervals by their length
    sortedData = data[:][(data[:,2]-data[:,1]).argsort()]
    piLengths = sortedData[:,2]-sortedData[:,1]

    # filter by the elbow method documented in this paper: Ville Satopa, Jeannie Albrecht, David Irwin, and Barath Raghavan,
    # 'Finding a “Kneedle” in a Haystack: Detecting Knee Points in System Behavior,' 31st International Conference on
    # Distributed Computing Systems Workshops, Minneapolis, MN, 2011, pp. 166-171, doi: 10.1109/ICDCSW.2011.20.
    if filterType == 'elbow' :
        if filterArg1 == None :
            filterArg1 = 1.0
        try :
            elbow = KneeLocator(range(sortedData.shape[0]), piLengths, S=filterArg1, curve='convex', direction='increasing')
            return applyCutoff(sortedData, elbow.elbow)
        except :
            # no elbow found
            return sortedData

    # filter out all PIs with lengths below a fixed value
    elif filterType == 'fixed' :
        if filterArg1 == None :
            sys.exit("The 'fixed' filter requires a cutoff length.  Aborting.")
        cutoff = bisect.bisect_left(piLengths, filterArg1)
        return applyCutoff(sortedData, cutoff)

    # filter by a scaled cutoff distance
    elif filterType == 'scaled' :
        if filterArg1 == None :
            sys.exit("The 'scaled' filter requires a cutoff length.  Aborting.")
        if filterArg2 == None :
            filterArg2 = .01
        filteredData = []
        for pi in sortedData :
            # keep if death > (birth + delta + (birth * scale))
            if pi[2] > (pi[1] + filterArg1 + (pi[1] * filterArg2)) :
                filteredData.append(pi)
        return np.array(filteredData)
        
    # filter out PIs with length below the max(0, mean - (n * stdev))
    elif filterType == 'mean' :
        if filterArg1 == None :
            filterArg1 = 1.0
        cutoff = bisect.bisect_left(piLengths, np.mean(piLengths) - (np.std(piLengths) * float(filterArg1)))
        return applyCutoff(sortedData, cutoff)

    # filter by a percentage of the longest PIs 
    elif filterType == '%longest' :
        # assume that the user accidently gave the percentage as 25 instead of .25
        if filterArg1 == None :
            sys.exit("The '%longest' filter requires a cutoff percentage.  Aborting.")
        if filterArg1 > 1 :
            filterArg1 = filterArg1 / 100.0
        cutoff = round(float(sortedData.shape[0]) * filterArg1)
        return applyCutoff(sortedData, cutoff)

    # filter by returning the longest N PIs
    elif filterType == 'longestN' :
        if filterArg1 == None :
            sys.exit("The 'longestN' filter requires a number to PIs to return.  Aborting.")
        return applyCutoff(sortedData, round(filterArg1))
        
    else :
        sys.exit('Abort: Filter type ({}) not supported'.format(filterType))

## separate the persistence intervals by dimension into a list (each list element contains the persistence intervals at some dimension)
def _pisByDim(data) :
    maxDim = int(np.max(data[:, 0]))
    filteredData = []
    for i in range(maxDim+1) :
        pisAtDim = data[np.where(data[:,0]== i)]
        if pisAtDim.shape[0] != 0 :
            filteredData.append(pisAtDim)
    return filteredData

def filterPIs(data, filter=('elbow', 1.0), byDim=True) :
    """
    Filter persistence intervals using one of several cutoff functions (generally computed on the sorted (decreasing) length of 
    each persisence interval).  The filtering can be performed dimension by dimension or as a common cutoff across all dimensions.
    If the cutoff function sets a cutoff number above the number of persistence intervals, the full set is returned (in this case,
    the set can either be all of the persistence intervals [byDim=False] or from one of the dimensions [byDim=True].

    Parameters
    ----------
    data : numpy.ndarray
        numpy [*, 3] matrix; each row is a persistence interval written as <dimension, birth, death>.

    filter : tuple=('str'{, float}*) setup filter type and a set of one or more parameters needed for the filter, default:
             ('elbow', 1.0); available filters with parameter defaults are:   

        ('elbow', sens=1.0) Use the knee finding algorithm documented in the paper: Ville Satopa, Jeannie Albrecht, David
            Irwin, and Barath Raghavan, 'Finding a “Kneedle” in a Haystack: Detecting Knee Points in System Behavior,' 31st
            Int. Conf. on Distributed Computing Systems, 166-171, 2011, doi: 10.1109/ICDCSW.2011.20.
            This is the default filtering function, the default sensitivity is 1.0 (as recommended in the paper). 

        ('fixed', len) Filter persistence intervals that are below a fixed length (len), no default

        ('scaled', delta, scale=0.01) (experimental) Filter persistence intervals that are scaled (scale) by birth to a fixed length
            (delta), no len default.  The the PI is kept if death > (birth + delta + (birth * scale)).

        ('mean', numStd=1.0) Filter persistence intervals by their mean + (numStd * standardDeviations).

        ('longestN', N) Filter persistence intervals by a count of the longest persistence intervals, no default.

        ('%longest', percent) Filter persistence intervals by a percentage of the longest persistence intervals, no default.


    byDim : boolean
        filter PIs by dimension (True) or filter all together (False) (default=True)

    Returns
    -------
    filtered persistence intervals


    Notes: 
    ------
    S the idea of filtering PIs is based on the common though that short PIs are mostly noise (although some studies do rely on
    those short PIs).  As a result, this filtering capability should be used carefully.  Several of the selection criteria are
    probably garbage (especially the latter set in the function); I implemented them because they occurred to me at some time in
    my studies with persistent homology.  I now believe that the most useful of these are the 'elbow', 'mean', and 'atLeast'
    filters.  I setting 'elbow' as the default cutoff as it should find the sharpest elbow in the sorted distances and remove
    (filter out) those below the elbow.  I am somewhat conflicted about this, as the 'atLeast' (remove all PIs with a length less
    than a specified distance....which should be short) is probably also a very powerful filtering mechanism; but I just cannot
    decide what distance should be the default and technically, I suspect it is probably an artifact of the data.  That is why, I
    also suspect that the mean plus various standard deviations ('mean') is also a strong candidate for use.
    """

    if byDim == True :
        filteredData = []
        for i in _pisByDim(data) :
            filteredData.append(_filterPIs(i, filter))
        return np.vstack((filteredData))
    else :
        return _filterPIs(data, filter)
