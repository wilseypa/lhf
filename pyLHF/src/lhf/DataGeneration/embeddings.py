#!/usr/bin/python3

## Embedding a d-dimensional object into an m-dimensional spaces (m >= d)

import sys
import numpy as np
from scipy.stats import ortho_group

def embedding(ptCld, targetDim = 3, offset = None, debug = False):
	"""
	Generates a point cloud in m dimensions from a point cloud in n dimensions. Naively embeds the point cloud in the higher
	dimensional space, then performs an orthogonal transformation using a random orthogonal matrix. 
	An orthogonal matrix will preserves distances and angles while rotating and reflecting the point cloud.

	Parameters
	----------
	ptCld : dict 
		ptCld["points"] : The point cloud
	targetDim : int 
	    The dimension of the space we will embed the point cloud in (default 3)
	offset : list (of length targetDim)
	    Vector to translate the point cloud by (default (0, 0, ..., 0))
	debug : boolean 
	    Not currently used (default False)

	Returns
	-------
	rotNewPtCld : A point cloud in dimension targetDim
	"""

	#Get the number of points and the dimension of the original point cloud
	n, d = ptCld['points'].shape

	if targetDim < d:
		print('Aborting: Can only embed a point cloud into the same or higher dimensional space.  Source Dim: {}, Target Dim: {}'.format(d,targetDim))
		sys.exit(-1)

	if offset is None:
		#By default, don't translate the point cloud
		offset = [0 for i in range(targetDim)]

	#Embed the point cloud in the target dimension by setting all additional axes to be 0
	newPtCld = np.zeros((n, targetDim))
	newPtCld[:,:d] = ptCld['points']

	#Translate the point cloud by the offset
	newPtCld = np.add(newPtCld, offset)

	#Generate a random orthogonal marix to rotate the point cloud
	rot = ortho_group.rvs(dim = targetDim)
	#Transform the point cloud using matrix multiplication
	rotNewPtCld = np.dot(newPtCld, rot)

	#Update the values in the dictionary to reflect the embedded point cloud
	#Pairwise distances are preserved, so we don't need to update any of the distance keys in the dictionary
	ptCld['ptCldObject'] = 'embedded ' + ptCld['ptCldObject']
	ptCld['dimension'] = targetDim
	ptCld['origin'] = np.zeros(targetDim)
	ptCld['geometricCenter'] = rotNewPtCld.mean(axis=0)
	ptCld['points'] = rotNewPtCld

	return ptCld
