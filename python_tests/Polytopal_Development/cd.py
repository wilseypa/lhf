
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import tadasets
from mpl_toolkits.mplot3d import Axes3D


#Schlegel diagram and Steriographic Projections for curvature removal as we move to lower dimensions
dimension = 3
points = tadasets.dsphere(n=100, d=dimension-1, r=1, noise=0)
points=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1],[0,0,0])

#Grid Precision as percentage of maximum diagonal/distance.
precision = .5
maximumdistance = 0;
distancematrix = []	
minimumdistance = 99999;
def computedistancematrix():
	global minimumdistance
	global maximumdistance
	for x in range(0,len(points)):
		distance = []
		for y in range(x,len(points)):
			di = math.dist(points[x],points[y])
			distance.append(di)
			if(di>maximumdistance):
				maximumdistance = di
			if(di<minimumdistance and di != 0):
				minimumdistance = di
		distancematrix.append(distance)
	return distancematrix

def distanceindex(x,y):
	distancematrix = computedistancematrix()
	if(x>y):
		return distancematrix[y][x-y]
	else:
		return distancematrix[x][y-x]
		
def gridedpoints(points,precison):
	gridprecision = precision*minimumdistance
	adjustedpointongrid = (points // gridprecision)*gridprecision+(gridprecision/2)
	return adjustedpointongrid
gridedpoints = gridedpoints(points,precision)	
computedistancematrix()	
tri = Delaunay(points)
print(len(tri.simplices))
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))
 
def union(lst1,lst2):
	return list(set().union(lst1 , lst2))

def neighbours(polytop,simplices):
	polytop = list(polytop)
	adjacency = []
	for y in simplices:
		y = list(y)
		if(polytop != y):
			if(len(intersection(polytop,y))==dimension):
				adjacency.append(y)
	return adjacency
#xt = int(math.sqrt(len(simplices)))+1
#fig, axs = plt.subplots(xt,xt)

def mergeneighbors(polytop,simplices):
	neighbors = neighbours(polytop,simplices)
	for x in neighbors:
		y = union(x,polytop)
		hull = ConvexHull([points[i] for i in y])
		hullboundary = {}
		for sublist in hull.simplices:
			for item in sublist:
				if item in hullboundary.keys():
					hullboundary[item] += 1
				else:
					hullboundary[item] = 1
		if(len(list(hullboundary.keys())) == len(y)):
			polytop = union(polytop,x)
	return polytop


#for in adjacency
def optimalconvexization(simplices):
	convexparts = []
	maxdist = []
	averagedist = []
	sizecd = []
	for i in range(0,len(simplices)):
		#axs[int(i/xt)][i%xt].scatter(points[:,0],points[:,1])
		x= list(simplices[i])
		distance =0
		maxval = 0
		while True:
			x1 = list(mergeneighbors(x,simplices))
			if(x1==x):
				break
			x = x1
		if x not in convexparts:
			maxval=0
			distance = 0
			convexparts.append(x)
			hull = ConvexHull([points[i] for i in x])
			for p in hull.simplices:
				distance = distance + math.dist(points[x[p[0]]],points[x[p[1]]])
				if(maxval < math.dist(points[x[p[0]]],points[x[p[1]]])):
					maxval = math.dist(points[x[p[0]]],points[x[p[1]]])
			maxdist.append(maxval)
			averagedist.append(distance)
			sizecd.append(len(x))
		k=set()
		for tp in convexparts:
			for lm in tp:
				k.add(lm)
		if(len(k)==len(points)):
			break;
	df = pd.DataFrame(list(zip(convexparts,maxdist,averagedist,sizecd)),columns =['convexpart', 'maxdist','averagedist','sizecd'])
	df = df.sort_values(by = 'sizecd',ascending = False)
	df = df.sort_values(by = 'maxdist',ascending = True)
	return df


valid = []
def iterativeconvexization(simplices):
	convexpartsunion = []
	simplices = list(simplices)
	remainingsimplices = simplices
	premaining = 0
	while(True):
		df = optimalconvexization(remainingsimplices)
		i=0
		pointsaddressed = []
		for x in df['convexpart'].tolist():
			check = 1
			for t in convexpartsunion:
				if(len(intersection(t,x)) > dimension):
					check = 0
			if(check==1):
				valid.append(i)
				hull = ConvexHull([points[i] for i in x])
				#for p in hull.simplices:
				#	axs[int(i/xt)][i%xt].plot([points[x[p[0]]][0],points[x[p[1]]][0]],[points[x[p[0]]][1],points[x[p[1]]][1]])
				pointsaddressed = union(pointsaddressed,x)
				convexpartsunion.append(x)
			i = i+1
		remainingsimplices = []
		remaining = 0
		for i in simplices:
			present =0
			for x in convexpartsunion:
				if (len(intersection(x,i))>=len(i)):
					present = 1
					break
			if(present!=1):
				remaining = remaining+1
				remainingsimplices.append(i)	
		print(remaining)
		if(remaining==premaining):
			if(remaining !=0):
				for r in remainingsimplices:
					convexpartsunion.append(r)
			break
		premaining = remaining
	return convexpartsunion

convexpartsunion = iterativeconvexization(tri.simplices)


X=np.array(points)
fig = plt.figure(figsize=(10,7))
ax=Axes3D(fig)
#ax.plot([i,i],[j,j],[k,h],color = 'g')
#plt.show()
X1 = np.transpose(X)[0]
Y1 = np.transpose(X)[1]
Z1 = np.transpose(X)[2]


for x in convexpartsunion:
	print(len(x))
	hull = ConvexHull([points[i] for i in x])
	for p in hull.simplices:
		for t in range(0,len(p)):
			for e in range(t+1,len(p)): 
				#plt.plot([points[x[p[0]]][0],points[x[p[1]]][0]],[points[x[p[0]]][1],points[x[p[1]]][1]])
				ax.plot([X1[x[p[t]]],X1[x[p[e]]]], [Y1[x[p[t]]],Y1[x[p[e]]]], [Z1[x[p[t]]],Z1[x[p[e]]]])
				 
#plt.scatter(points[:,0],points[:,1],s=10,color='black')
ax.scatter(np.transpose(X)[0], np.transpose(X)[1], np.transpose(X)[2])



convexp = iterativeconvexization(hull.simplices)

plt.show()
'''
for x in convexpartsunion:
	print(x)
	hull = ConvexHull([points[i] for i in x])
	for p in hull.simplices:
		print(p)

'''
