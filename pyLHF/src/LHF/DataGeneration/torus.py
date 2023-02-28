#!/usr/bin/python3

## This is a simple generator to build torii for testing TDA tools/techniques. 

import sys
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
import math
import random
from ..utilities import rejection_sample, run_trials
from .dSphere import dSphere_uniformRandom
import datetime

def torus_noisy_surface_uniformRandom(dimension=2, numPoints=1000, r1=1.0, r2=2.0, trials=16, selection='max') :
    """
    Generate a point cloud for the boundary of a torus populated with randomly distributed points.

    This function generates a torus with intrinsic dimension {dimension} embedded in R^{2 * dimension}. For example, a
    2-torus will be embedded in R^4, and this is known as a Clifford Torus: see https://en.wikipedia.org/wiki/Clifford_torus.
    This can be generalized to higher dimensions, as described in the wikipedia pages, and each such torus is "flat" in the sense
    that is has zero Gaussian curvature: see https://en.wikipedia.org/wiki/Torus#Flat_torus.

    Even though this torus is embedded in R^{2 * dimension}, its maximum homology group should only be H_{dimension}. The Betti 
    numbers for this torus should follow the binomial coefficients; for a 2-torus in R^4, these should be 1 in H0, 2 in H1, and 1 in H2.
    Forr a 3-torus, this will be 1 in H0, 3 in H1, 3 in H2, and 1 in H3.

    While a d-dimensional torus is typically embedded in R^{dimension+1}, the Clifford Torus is embedded in R^{2*dimension}. It
    is generated from the Cartesian Product of dimension S^1 circles. A point can be sampled from an S^1 
    circle by uniformly choosing an angle between 0 and 2pi with radius r in polar coordinates. A point on the torus can then be sampled as
    ((circle 1)_x, (circle 1)_y, ..., (circle d)_x, (circle d)_y))

    To generate a thick boundary around the surface of the torus, we sample points uniformly from an annulus instead of a circle.

    When dimension = 1, this generates a circle in R^2
    WHen dimension = 2, this generates a T^2 torus in R^4. This is the same torus that we typically think of when we refer to a torus, but
    embedded in R^4 instead of R^3

    Parameters
    ----------
    dimension : int 
        The number of dimensions for the desired torus (default 2)
    numPoints : int 
        The number of points desired in the output point cloud (default 1000)
    r1 : either a float (default) or list of floats of length dimension
        [float]: The ith entry is the inner radius of the ith circle
        float: The radius of the point cloud (default 1.0). Every circle will have the same radius
    r2 : either a float (default) or list of floats of length dimension
        [float]: The ith entry is the outer radius of the ith circle
        float: The radius of the point cloud (default 2.0). Every circle will have the same radius

    trials: int 
        The number of trials used to optimize the distribution of points in the dSphere (default: 16)
    selection: str (min, max, mean, var)
        The selection criteria to select between two candidate point clouds:
            'min' : maximize the minimum nearest neighbor distance of the points
            'max' (default): minmimize the maximum nearest neighbor distance of the points
            'mean' : maximize the mean of the nearest neighbor distances.
            'stdev' : minimize the standard deviation of the nearest neighbors

    Returns
    -------
    torus: dict

        torus['ptCldObject'] = 'Clifford Torus'
        torus['origin'] = (0, 0, ...) : 
        torus['radius'] = (r1, r2) : the radius of the point cloud boundary
        torus['points'] = randomPoints : the point cloud

    """

    #### --------------------------------------------------------------------------------
    #### helper functions
    #### --------------------------------------------------------------------------------
    
    # Sample a point uniformly on the inside of an annulus with inner radius r1 and outer radius r2 
    def uniform_circle_point(r1, r2):
        if r1 > r2:
            print("Can't generate a circle with inner radius greater than the outer radius")
            sys.exit(-1)
        theta = 2 * random.random() * math.pi
        r = math.sqrt(random.random() * (r2**2 - r1**2) + r1**2)
        return r * math.cos(theta), r * math.sin(theta)

    # Sample a point uniformly on the boundary of a torus with inner radii r1 and outer radii r2 in dimension d
    def uniform_torus_point(d, r1, r2):
        #Independently for each i, sample coordinates x_{2i-1} and x_{2i} from the boundary of a circle
        return [coord for i in range(d) for coord in uniform_circle_point(r1[i], r2[i])]

    # Generate n points on the boundary of a torus with radius r in dimension d
    def get_torus(dimension, numPoints, r1, r2):
        if dimension <= 0 :
            print("No 0-torii exist. Aborting.`")
            sys.exit(-1)

        if type(r1) == list:
            if len(r1) != dimension:
                print("A radius must be specified for each dimension of the torus")
                sys.exit(-1)
        else:
            #Specify that the radius of each circle is equal
            r1 = [r1 for _ in range(dimension)]

        if type(r2) == list:
            if len(r2) != dimension:
                print("A radius must be specified for each dimension of the torus")
                sys.exit(-1)
        else:
            #Specify that the radius of each circle is equal
            r2 = [r2 for _ in range(dimension)]

        randomPoints = np.array([uniform_torus_point(dimension, r1, r2) for _ in range(numPoints)])
        
        return randomPoints
    
    #### --------------------------------------------------------------------------------

    sample_func = lambda: get_torus(dimension, numPoints, r1, r2)
    points = run_trials(sample_func, trials, selection)

    torus = dict()
    torus['ptCldObject'] = 'Clifford Torus'
    torus['date'] = datetime.datetime.now()
    torus['origin'] = np.zeros(2*dimension)
    torus['radius'] = (r1, r2)
    torus['points'] = points

    return torus

def torus_boundary_uniformRandom(dimension=2, numPoints=1000, r=1.0, trials=16, selection='max') :
    """
    Generate a point cloud for the boundary of a torus populated with uniformly distributed points.
    
    The parameters are the same as the torus_noisy_surface_uniformRandom function, but both r1 and r2 are set to r to 
    guarantee that the points are selected uniformly randomly.
    """
    return torus_noisy_surface_uniformRandom(dimension, numPoints, r, r, trials, selection)

def box_genus_with_boundary_uniformRandom(numPoints=1000, genus=1, trials=16, selection='max'):
    '''
    Generate a point cloud for the boundary of a surface with genus populated with uniformly distributed points. The
    surface is generated by sampling from the sides of a rectangular prism, and then inserting tunnels for each genus.
    The result is a box-like looking surface which is homeomorphic to the torus with corresponding genus.

    The points are in R3. The box has dimensions (2*genus + 1) x 3 x 1, and each tunnel has dimensions 1 x 1 x 1.
    A horizontal cross section of the double torus looks like the shape below.
    _________________
    |   __     __   |
    |  |__|   |__|  |
    |_______________|

    Parameters
    ----------
    numPoints: int 
        The number of points desired in the output point cloud (default 1000)
    genus: int
        The desired genus of the surface (default: 1)
    trials: int 
        The number of trials used to optimize the distribution of points in the dSphere (default: 16)
    selection: str (min, max, mean, var)
        The selection criteria to select between two candidate point clouds:
            'min' : maximize the minimum nearest neighbor distance of the points
            'max' (default): minmimize the maximum nearest neighbor distance of the points
            'mean' : maximize the mean of the nearest neighbor distances.
            'stdev' : minimize the standard deviation of the nearest neighbors

    Returns
    ----------
    torus: dict

        torus['ptCldObject'] = 'Box Torus'
        torus['origin'] = (0, 0, 0) : 
        torus['genus'] = genus
        torus['points'] = randomPoints : the point cloud
    '''

    def sampleTorus(n, g):
        # Sample points from the interior of the box 
        pts = np.column_stack((np.random.uniform(0, 2*g+1, n), np.random.uniform(0, 3, n), np.random.uniform(0, 1, n)))
        rand = np.random.randint(0, 20*g+14, n)
        for i in range(n):
            # Weight by surface area of sides to ensure uniform random sampling
            if rand[i] < 12*g + 6:
                # Move point to top or bottom
                pts[i][2] = random.randint(0,1)
            elif rand[i] < 16*g + 8:
                # Move point to left or right side
                pts[i][1] = random.randint(0, 1) * 3
            elif rand[i] < 16*g + 14:
                # Move point to front or back side
                pts[i][0] = random.randint(0, 1) * (2*g + 1)
            elif rand[i] < 18*g + 14:
                # Add tunnel front and back sides
                pts[i][0] = random.randint(1, 2*g)
                pts[i][1] = random.uniform(1, 2)
            else:
                # Add tunnel left and right sides
                pts[i][0] = random.uniform(0, 1) + random.randint(0, g-1)*2 + 1
                pts[i][1] = random.randint(1, 2)
        return pts

    def acceptTorus(points):
        # Remove points inside the tunnels
        return (points[:,0] % 2 <= 1) | (points[:,1] <= 1) | (points[:,1] >= 2)

    sample_func = lambda: rejection_sample(numPoints, lambda numPoints: sampleTorus(numPoints, genus), acceptTorus)
    points = run_trials(sample_func, trials, selection)

    torus = dict()
    torus['ptCldObject'] = 'Box Torus'
    torus['genus'] = genus
    torus['date'] = datetime.datetime.now()
    torus['origin'] = np.zeros(3)
    torus['points'] = points

    return torus

def box_genus_interior_uniformRandom(numPoints=1000, genus=1, trials=16, selection='max'):
    '''This function is exactly analogous to box_genus_with_boundary_uniformRandom, but it returns points sampled
    uniformly randomly from the interior of a surface with genus g instead of the boundary
    '''
    def sampleSolidTorus(n, g):
        return np.column_stack((np.random.uniform(0, 2*g+1, n), np.random.uniform(0, 3, n), np.random.uniform(0, 1, n)))

    def acceptSolidTorus(points):
        return (points[:,0] % 2 < 1) | (points[:,1] < 1) | (points[:,1] > 2)

    sample_func = lambda: rejection_sample(numPoints, lambda numPoints: sampleSolidTorus(numPoints, genus), acceptSolidTorus)
    points = run_trials(sample_func, trials, selection)

    torus = dict()
    torus['ptCldObject'] = 'Box Torus Interior'
    torus['genus'] = genus
    torus['date'] = datetime.datetime.now()
    torus['origin'] = np.zeros(3)
    torus['points'] = points

    return torus

def spherical_genus_boundary_random(numPoints = 1000, genus = 3, radius = .3, trials = 16, selection = 'max'):
    """
    Generate a point cloud for the boundary of a torus with genus 1, 3, or 5. These points will NOT be distributed 
    uniformly randomly along the boundary. 

    This function first samples points uniformly random from the surface of a sphere with radius 1. The function then
    adds tunnels to add genus. For genus 1, a single tunnel is added by sampling points from a cylinder.
    This resembles the standard torus point cloud. For genus 3, a second cross-tunnel is added
    (reference Figure 3 of "Estimating Betti Numbers Using Deep Learning"). For genus 5, a third
    cross-tunnel is added (reference https://mathoverflow.net/questions/98925/building-a-genus-n-torus-from-cubes)

    Parameters
    ----------
    numPoints: int 
        The number of points desired in the output point cloud (default 1000)
    genus: (1, 3, or 5)
        The desired genus of the surface (default 3)
    radius:
        The radius of each tunnel (default .3)
    trials: int 
        The number of trials used to optimize the distribution of points in the dSphere (default: 16)
    selection: str (min, max, mean, var)
        The selection criteria to select between two candidate point clouds:
            'min' : maximize the minimum nearest neighbor distance of the points
            'max' (default): minmimize the maximum nearest neighbor distance of the points
            'mean' : maximize the mean of the nearest neighbor distances.
            'stdev' : minimize the standard deviation of the nearest neighbors

    Returns
    ----------
    torus: dict

        torus['ptCldObject'] = 'Spherical Torus'
        torus['origin'] = (0, 0, 0)
        torus['genus'] = genus
        torus['points'] = randomPoints : the point cloud
    """
    def sampleCylinder(n, r, h):
        """Sample n points from a cylinder of radius r and heigh h (z-values stretching from -h/2 to h/2)"""
        angles = np.random.uniform(0, 2*math.pi, n)
        return np.column_stack((r*np.cos(angles), r*np.sin(angles), np.random.uniform(-h/2, h/2, n)))

    def sampleSphericalTorus(n, g, r):
        if g not in [1, 3, 5]:
            raise "Invalid genus: only 1, 3, or 5 supported"
        cylinderN = n // 8
        sphereN = n - cylinderN*((g+1)//2)
        sphere = dSphere_uniformRandom(dimension=3, numPoints=sphereN, r1=1, r2=1, trials=1)['points']
        cylinder = sampleCylinder(cylinderN, r, 2)
        # Add the first tunnel in the z direction
        points = np.concatenate((sphere, cylinder))
        if g == 3 or g == 5:
            # Add a cross tunnel in the x direction
            cylinder = sampleCylinder(cylinderN, r, 2)[:,[2,0,1]]
            points = np.concatenate((points, cylinder))
        if g == 5:
            # Add a cross tunnel in the y direction
            cylinder = sampleCylinder(cylinderN, r, 2)[:,[0,2,1]]
            points = np.concatenate((points, cylinder))
        np.random.shuffle(points)
        return points

    def acceptSphericalTorus(points, g, r):
        mask = (points[:,0]**2 + points[:,1]**2) >= r ** 2
        if g > 1:
            mask &= (points[:,1]**2 + points[:,2]**2) >= r ** 2
        if g == 5:
            mask &= (points[:,0]**2 + points[:,2]**2) >= r ** 2
        return mask
        
    sample_func = lambda: rejection_sample(numPoints, lambda numPoints: sampleSphericalTorus(numPoints, genus, radius), 
                                        lambda points: acceptSphericalTorus(points, genus, radius))
    points = run_trials(sample_func, trials, selection)

    torus = dict()
    torus['ptCldObject'] = 'Spherical Torus'
    torus['genus'] = genus
    torus['date'] = datetime.datetime.now()
    torus['origin'] = np.zeros(3)
    torus['points'] = points

    return torus