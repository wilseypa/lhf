# Ref. https://github.com/outlace/OpenTDA
# A set of functions for building a Vietoris-Rips complex from point cloud data and
# then collapse it into a folded complex.

import numpy as np
import matplotlib.pyplot as plt

# Euclidian distance function.
def euclidianDist(x, y):
    return np.linalg.norm(x - y)

# Build neighorbood graph G which is the 1-skeleton.
def build1Skeleton(dataSet, epsilon = 5, dist = euclidianDist):   # dataSet is a numpy array
    # Count the number of rows (points) in the data set. Create an array, called nodes, that
    # has all the points of the data set represented as [0, 1, 2, ..., n], where n is the 
    # number of points in the data set.
    nodes = [points for points in range(dataSet.shape[0])]
    edges = []     # Create an empty edge array.
    # Create an empty weight array that will store the weight or length of each edge.
    weights = []
    for i in range(dataSet.shape[0]):   # Iterate through each data point
        for j in range(dataSet.shape[0] - i):   # Inner loop to calculate pairwise point distances between data points
            # Iterate over the combination of each pair of points.
            a = dataSet[i]
            b = dataSet[i+j]
            if (i != j+i):   # Barring the combination of a point with itself
                distance = dist(a, b)
                if distance <= epsilon:
                    edges.append({i, i+j})   # Add an edge if the distance between points <= epsilon.
                    weights.append([len(edges)-1, dist])   # Store the weight of the edge preceded by an index.
                    
    return nodes, edges, weights


'''
We shall now expand the 1-skeleton into a k-complex according to the incremental algorithm
given by zomorodian-10. In order to do that, we assume a total (arbitrary) ordering on the
nodes in the neighborhood graph G. We begin with a simple utility function that finds all
neighbors of a node u within G that precede it in the given ordering.
'''
def lowerNbrs(nodes, edges, u):
    return {v for v in nodes if (u > v and {u,v} in edges)}