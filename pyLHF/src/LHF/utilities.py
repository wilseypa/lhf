#!/usr/bin/python3

## a collection of support utilities for our tda explorations.

from random import sample
import sys
import os
import numpy as np
from sklearn.random_projection import GaussianRandomProjection 
from sklearn.random_projection import SparseRandomProjection 
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from scipy.spatial import distance_matrix
import pylab
from mpl_toolkits.mplot3d import Axes3D


def cliProgressBar(iteration, totalIterations):
	"""
	A quick progress bar for some indicator while running
	
	Parameters
    ----------
    iteration : int
        Current iteration number
    totalIterations: int
        Total number of iterations

    """
    
	percent = ("{0:.1f}").format(100 * (iteration / float(totalIterations)))
	filled = int(100 * iteration // totalIterations)
	bar = '█' * filled + '-' * (100 - filled)
	print(f'\rProgress: |{bar}| {percent}% Complete', end='' )
	
	if iteration == totalIterations:
		print('\n')
	
	return
	

def embedPointCloud(ptCld, dim, rot='yes', fillVals=(0.0)) :
    """
    Parameters
    ----------
    ptCld : numpy.array(n, d)
        point cloud in dimension d
    dim: int
        dimension of the point cloud to embed the original point cloud within
    rotate : str
        'yes' (default), perform a random rotation of the high dimensional space (this will rotate the embedded point cloud in the
        higher dimensional space (the d-dimensional shape will be there, but its orientation will be randomly positioned);
        !='yes', no rotation is performed on the final point cloud 
    fillVals : list(int)
        list of fill values to use for the additional dimensions (default: (0.0))

    Return
    ------
    numpy.array(n, dim)

    """
    size = ptCld.shape

    if size[1] < dim :
        print ('Embedding dimension (' + str(dim) + ') is less than the dimension of the point cloud (' + str(size[1]) + ').  Aborting....')
        sys.exit(-1)

    if len(fillVals) == 1 :
        fillVals = np.full((size[0], dim - size[1]), fillVals)

    newPC = np.concatenate(ptCld, fillVals)
    if rot == 'yes' :
        rotate, _ = np.linalg.qr(np.random.random((dim, dim)))
        return np.dot(newPC, rotate)
    else :
        return newPC
    
    
def projectDataTo3D(ptCld, projectionType='gaussian') :
    """ 
    Project high dimensional point cloud to a 3D point cloud.

    Perform random projection to achieve data reduction of a high dimensional point cloud to a point cloud in 3D.  The function
    uses either DB Friendly (SparseRandomProjection, with modifications from Ping-06-KDD) or JL (GaussianRandomProjection) methods
    for the data reduction step.

    Parameters
    ----------
    ptCld : numpy.array(n, d)
        Cartisian coordinates in R^d to project 
    projectionType : str
        type of projection to perform: 'gaussian' (default), following the JL lemma; or 'sparse', DB Friendly

    Returns
    -------
    numpy.array(n, 3)
        Cartisian coordinates of reduced data in R^3

    """
    if projectionType == 'gaussian' :
        projector = GaussianRandomProjection(n_components=3)
    elif projectionType == 'sparse' :
        projector = SparseRandomProjection(n_components=3)
    else :
        print ('testingUtilities.projectDataTo3D: illegal projection projection type requsted; only gaussian or sparse supported')
        sys.exit(1)
    return projector.fit_transform(ptCldList)



## not sure this method is necessary, but i'm leaving it in here since its here....we may want to remove it in the figure.
def normalize(ptCld, method='MinMaxScalar') :
    """
    Normalize a point cloud.

    Normalize a point cloud using the sklearn functions MinMaxScalar[-1.0,1.0] or StandardScalar.

    Parameters
    ----------
    ptCld : numpy.array(n, d) 
        Point cloud in R^d to normalize
    method : string
        The scaling method from sklearn to use: 'MinMaxScalar' (the default) or StandardScalar

    Returns
    -------
    numpy.array(n, d)
        Normalized point cloud in R^d


    """
    if method == 'StandardScalar' :
        standardizer = StandardScaler().fit(ptCld)
        return standardizer.transform(ptCld)
    elif method == 'MinMaxScalar' :
        scaler = MinMaxScaler(feature_range=(-1.0,1.0))
        return scaler.fit_transform(ptCld)
    else :
        print ('testingUtilities.normalize: illegal normalization type requsted; only StandardScalar or MinMaxScalar supported')
        sys.exit(1)
        return ptCld



def plot3DPointCloud(ptCld, title='', setAxis=True) :
    """
    Generate a 3D plot of a point cloud.

    Parameters
    ----------
    ptCld : numpy.array(n, 3) 
        Cartisian coordinates in R^3 to plot
    title : string
        String to place in figure title, if blank no title generated
    setAxis : boolean
		If true, scales all axis to min/max of PC data; otherwise default
    """
    figure = pylab.figure()
    fig = figure.add_subplot(111, projection='3d')
    if title != '' :
        pylab.title(title)
    #fig.set(aspect='equal')
	
	#NOTE: Currently set this behind a flag but for PCs not normalized the stretch needs to be updated to take the 
	#       maximum difference of any dimension as the bounding box (this will be tighter and keep aspect ratio)
    if(setAxis):
        max = np.max(ptCld)
        min = np.min(ptCld)
        # i prefer the same boundary limits on the axes
        if (min < 0) :
            if -min > max :
                max = -min
            else :
                min = -max
        fig.set_xlim3d(min, max)
        fig.set_ylim3d(min, max)
        fig.set_zlim3d(min, max)
        #fig.set_axis_off()
        fig.set_xticks([min, 0, max])
        fig.set_yticks([min, 0, max])
        fig.set_zticks([min, 0, max])

    fig.scatter(ptCld[:,0], ptCld[:,1], ptCld[:,2], marker='.', s=3)
    #pylab.show()
    return figure
    
# this next function is planned for removal

## read a list of csv files containing point cloud data.  each file is read into a numpy matrix and stored in a standard
## python list.  the rows of the csv files each represent a unique n-dimensional point in the n-space.  each csv file
## can contain points for a a unique n-space. 
def inputDataFiles(fileNameList) :
    ptCld = []

    for i in fileNameList :
        # add .csv suffix if missing
        if not(i.endswith('.csv')) :
            i += '.csv'
        ptCld.append(np.loadtxt(i, delimiter=',', comments='#', dtype=np.float))

    return ptCld




# Referenced from https://stackoverflow.com/a/50107963
def rejection_sample(n, sample_func, accept_func):
    '''Perform rejection sampling to obtain n valid points.
    
    Parameters
    ----------
    n: int
        The number of points to sample
    sample_func: f(n) -> points
        A function with one parameter, n, which returns n points randomly
        sampled from a surface
    accept_func: f(points) -> mask
        A function with one parameter, points, which is an array of points,
        which returns a boolean array (true if the point is valid and false otherwise) 

    Return
    ------
    points: numpy.array(n, dim)
    '''
    oversample = 5 * n  # Arbitrary oversample by 5 to limit the number of calls to sample_func
    points = sample_func(oversample)
    mask = accept_func(points)
    reject, = np.where(~mask)
    while oversample - reject.size < n:
        fill = sample_func(reject.size)
        mask = accept_func(fill)
        points[reject[mask]] = fill[mask]
        reject = reject[~mask]
    return np.delete(points, reject, axis = 0)[:n]

def run_trials(sample_func, trials, selection):
    '''
    Sample a point cloud multiple times and return the best sample according to a selection criteria

    Parameters
    ----------
    sample_func: f() -> np.array(numPoints, dim)
        Returns points sampled from a surface
    trials: int 
        The number of trials used to optimize the distribution of points
    selection: str (min, max, mean, var)
        The selection criteria to select between two candidate point clouds:
            'min' : maximize the minimum nearest neighbor distance of the points
            'max' (default): minmimize the maximum nearest neighbor distance of the points
            'mean' : maximize the mean of the nearest neighbor distances.
            'stdev' : minimize the standard deviation of the nearest neighbors'''
    points = sample_func()
    pointsNN = np.sort(np.sort(distance_matrix(points, points))[:,1])

    for i in range(trials) :
        pointsTrial = sample_func()
        pointsTrialNN = np.sort(np.sort(distance_matrix(pointsTrial, pointsTrial))[:,1])
        selectTrial = False

        if selection == 'mean' and np.mean(pointsTrialNN) > np.mean(pointsNN) : selectTrial = True
        if selection == 'stdev' and np.std(pointsTrialNN) < np.std(pointsNN) : selectTrial = True
        if selection == 'min' and pointsTrialNN[0] > pointsNN[0] : selectTrial = True
        if selection == 'max' and pointsTrialNN[-1] < pointsNN[-1] : selectTrial = True
        if selectTrial == True :
            points = pointsTrial
            pointsNN = pointsTrialNN
 
    return points
