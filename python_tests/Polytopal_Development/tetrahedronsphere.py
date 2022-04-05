from sklearn.decomposition import PCA
import math
import numpy as np
import pandas as pd
from itertools import combinations
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.axes(projection='3d')


dim = 3
edgecount = (int)(dim*(dim+1)/2)
points = 500+edgecount+dim+1
pointperedge = (int)(points/edgecount)
poly = [i for i in range(0,dim+1)]
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
	

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
        
for dimension in range(dim,dim+1):
	points = []
	pnts = dimension+1
	cnt  =0
	origin = [0 for x in range(0,dimension)]
	for x in range(0,pnts):
		coords = []
		for y in range(0,dimension+1):
			if(y==cnt):
				coords.append(1)
			else:
				coords.append(0)
		points.append(coords)
		cnt = cnt+1
    
	newpoints = PCAprojection(points)
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
ax.set_zlim(-1,1)
plt.savefig("tetrasphere.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()
