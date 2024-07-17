#!/usr/bin/python3

## a collection of triangulation utilities for our tda explorations.

import sys
import os
import numpy as np
from scipy.spatial import Delaunay
from itertools import combinations
from scipy.sparse import csr_matrix

#define a dynamic object for listing
class Object(object):
	pass

def getTriangulation(data):
	"""
	Get the delaunay triangulation of the input data. Use the maximum edge weight for
		creating simplices, then order the simplices by weight for filtration.
		
    Parameters
    ----------
	data : numpy.array(n, d)
		point cloud for (d+1)-triangulation
		
		
	Return
    ------
    numpy.array(n, dim)
		Simplex array list for (d+1)-triangles with maximal weights
	"""

	a = Delaunay(data, qhull_options="QJ")

	Kf = []
		
	# Insert simplices and maximum edge weight
	for i in range(len(a.simplices)):

		##Set to the first distance
		maxDist = np.linalg.norm((data[a.simplices[i,0]])-(data[a.simplices[i,1]]))
		i_simp = [(a.simplices[i,0], a.simplices[i,1])] 

		#Check remaining simplices with combinations
		for it_simp in combinations(a.simplices[i], 2):
			cur = np.linalg.norm(data[it_simp[0]] - data[it_simp[1]])
			if cur > maxDist:
				maxDist = cur
				i_simp = it_simp

		t = Object()
		t.weight = maxDist

		t.simplex = a.simplices[i]
		t.simplex.sort()
		Kf.insert(0,t)

	# Sort d-simplices by weight
	Kf.sort(key=lambda x: x.weight)

	return Kf


def projectTriangulation(tri, dmat, d):
	"""
	Create a new triangulation by decomposing the original (tri) into dimension d.
		
    Parameters
    ----------
	tri : numpy.array(n, d)
		Simplex array list for (d+1)-triangles with maximal weights
		
	dmat : numpy.array(n, d)
		Distance matrix for the original point clouod
		
	d : int
		Dimension of d-triangles to generate
		
	Return
    ------
    numpy.array(n, dim)
		Simplex array list for (d)-triangles with maximal weights
	"""
	newTri = []
	insertedSets = set()
	
	for tri_i in range(len(tri)):
		
		for it_simp in combinations(tri[tri_i].simplex, d):
		
			if(it_simp in insertedSets):
				continue;
		
			maxDist = dmat[it_simp[0],it_simp[1]]
			i_simp = [(it_simp[0], it_simp[1])]
			
			for vert in combinations(it_simp, 2):
				
				cur = dmat[vert[0], vert[1]]
				if cur > maxDist:
					maxDist = cur
					i_simp = it_simp
					
			t = Object()
			t.weight = maxDist
			t.simplex = it_simp
			newTri.insert(0, t)
			insertedSets.add(it_simp)
	
	newTri.sort(key=lambda x: x.weight)
	
	return newTri


def triToSparse(tri, data):
	"""
	Create a sparse matrix for the triangulation data (for PH approximation)
		
    Parameters
    ----------
	tri : numpy.array(n, d)
		Simplex array list for (d+1)-triangles with maximal weights
		
	data : numpy.array(n, d)
		Original point cloud data
		
	Return
    ------
    numpy.array(n, dim)
		Sparse incidence matrix based on triangulation
	"""
	#create np sparse distance matrix
	sparseMat = np.zeros((len(data), len(data)))
	
	for i in range(len(tri)):
		for it_simp in combinations(tri[i].simplex, 2):
			sparseMat[it_simp[0], it_simp[1]] = tri[i].weight
			sparseMat[it_simp[1], it_simp[0]] = tri[i].weight
			#sparseMat[it_simp[0], it_simp[1]] = np.linalg.norm((data[it_simp[0]]) - (data[it_simp[1]]))
			#sparseMat[it_simp[1], it_simp[0]] = np.linalg.norm((data[it_simp[0]]) - (data[it_simp[1]]))
	
	sparseMat = csr_matrix(sparseMat)
	
	return sparseMat
