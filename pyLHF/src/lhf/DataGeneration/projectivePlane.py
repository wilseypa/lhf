#!/usr/bin/python3

## This is a simple generator to build real projective planes for testing TDA tools/techniques. 

import numpy as np
from .dSphere import dSphere_uniformRandom
from ..utilities import run_trials
import datetime

def real_projective_plane(numPoints = 20000, trials = 8, selection = 'min'):
    '''
    There is an embedding of RP2 in R4 defined by sampling points from a sphere and mapping (x, y, z) 
    # to (xy, xz, y^2 − z^2, 2yz). This embedding is defined in https://en.wikipedia.org/wiki/Real_projective_plane. 
    # However, this sampling is not guaranteed to be uniform.

    Parameters
    ----------
    numPoints: int 
        The number of points desired in the output point cloud (default 20000)
    trials: int 
        The number of trials used to optimize the distribution of points in the dSphere (default: 8)
    selection: str (min, mean, var)
        The selection criteria to select between two candidate point clouds:
            'min' : maximize the minimum nearest neighbor distance of the points
            'mean' : maximize the mean of the nearest neighbor distances.
            'stdev' : minimize the standard deviation of the nearest neighbors

    Returns
    -------
    plane: the point cloud
    '''

    def getPlane(numPoints):
        sphere = dSphere_uniformRandom(3, numPoints, r1=1.0, r2=1.0, trials=1)
        plane = np.zeros((numPoints, 4))
        plane[:,0] = sphere[:,0] * sphere[:,1]
        plane[:,1] = sphere[:,0] * sphere[:,2]
        plane[:,2] = sphere[:,1]**2 - sphere[:,2]**2
        plane[:,3] = sphere[:,1] * sphere[:,2]
        return plane

    sample_func = lambda: getPlane(numPoints)
    return run_trials(sample_func, trials, selection)
