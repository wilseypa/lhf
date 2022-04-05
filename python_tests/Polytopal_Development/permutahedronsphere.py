from sklearn.decomposition import PCA
import math
import numpy as np
import pandas as pd
from itertools import combinations
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import itertools

fig = plt.figure()
ax = plt.axes(projection='3d')


dim = 4
b = [i for i in range(0,dim)]

vertices = 24
edgecount = vertices*(vertices-1)/2
points = 1500+edgecount+dim+1
pointperedge = (int)(points/edgecount)
print(pointperedge)
poly = [i for i in range(0,24)]

def genPermutahedron(a, size):
    if size == 1:
        inputpoints.append(np.array(a))
        return
 
    for i in range(size):
        genPermutahedron(a, size-1)
        if size & 1:
            a[0], a[size-1] = a[size-1], a[0]
        else:
            a[i], a[size-1] = a[size-1], a[i]
            

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
        
def PCAprojection(pts):
	pca = PCA(n_components=len(pts[0])-1).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x

def generatesimplexfacets(poly,dimen):  #generate simplices of dimension dimen
	tp = list(combinations(poly, dimen))
	listreturn = []
	for x in tp:
		listreturn.append(list(x))
	return listreturn
	
for dimension in range(dim-1,dim):
	pnts = dimension+1
	cnt  =0
	origin = [0 for x in range(0,dimension)]
	inputpoints = []
	size = len(b)
	genPermutahedron(b,size)	
	pca_x = PCAprojection(inputpoints)
	inputpoints = pca_x
	inputpoints = np.array(inputpoints)
	newpoints = inputpoints
	points = set()
	edges = generatesimplexfacets(poly,2)
	for e in edges:
		for t in np.linspace(0, 1, pointperedge):
			points.add(tuple(t*(newpoints[e[0]])-newpoints[e[1]]*(t-1)))
	unitsphere = []
	for x in points:
		x = np.array(x)
		unitsphere.append(x/math.dist(origin,x))
axisEqual3D(ax)


zdata = [item[0] for item in unitsphere]
xdata = [item[1] for item in unitsphere]
ydata = [item[2] for item in unitsphere]
ax.scatter3D(xdata, ydata, zdata)
plt.grid(b=None)
plt.axis('off')
# Hide axes ticks
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1	,1)
plt.savefig("permutasphere.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()
