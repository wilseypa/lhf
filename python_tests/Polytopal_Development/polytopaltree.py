from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import tadasets
from mpl_toolkits.mplot3d import Axes3D
import polytope as pc
import pypoman.duality as polyfun
from sklearn.decomposition import PCA
from itertools import combinations
import heapq as hq
from functools import cmp_to_key
import csv
from random import randint
import copy

   
def minimumspanningtree(edges1):
	# using Kruskal's algorithm to find the cost of Minimum Spanning Tree
	res = 0
	pivots = []
	mstSize = 0
	ds = DisjointSet(datasize)
	for x in reversed(edges1):
		if ds.find(x.polytop[0]) != ds.find(x.polytop[1]):
			ds.union(x.polytop[0], x.polytop[1])
			res += x.weight
			mstSize +=1
			pivots.append(x)
			bettitableentry = [0,0,x.weight]
			bettieTable.append(bettitableentry)
			if(mstSize >= len(inputpoints)-1):
				for i in range(0,len(inputpoints)):
					if(ds.find(i) == i):
						bettitableentry = [0,0,'inf']
						bettieTable.append(bettitableentry)
				return pivots
	return pivots

def getDimEdges(dimension):
	return orderedpolytoparraylist[dim -dimension]

def getAllCofacets(e,dimension):    #will update this function a little
	returnlist = []
	for x in orderedpolytoparraylist[dim-dimension-1]:
		if(len(intersection(x.polytop,e.polytop)) == len(e.polytop)):
			#if(len(x.polytop) == dimension+2):
			returnlist.append(x)
	return returnlist

'''			else:
				minweight = max(x.weight,e.weight)
				otheredges = [y for y in x.polytop if y not in e.polytop]
				for check in otheredges:
					maxweight = e.weight
					valid = 1
					for f in e.polytop:
						if(incidenceindexreturn(check,f)!=1):
							valid = 0
							break
					if(valid==1):
						for f in e.polytop:		
							if(maxweight < distanceindex(check,f)):
								maxweight = distanceindex(check,f)
					if(minweight>maxweight):
						minweight = maxweight
				x.weight = minweight
				returnlist.append(x)
'''
			
def letter_cmp(a, b):
    if a.weight > b.weight:
        return -1
    elif a.weight == b.weight:
        if a.polytop > b.polytop:
            return 1
        else:
            return -1
    else:
        return 1


letter_cmp_key = cmp_to_key(letter_cmp)
	
def persistenceByDimension( edges, pivots, dimension):
	pivotcounter = len(pivots)-1
	edges.sort(key = letter_cmp_key)
	pivots.sort(key = letter_cmp_key)
	nextPivots = []	
	v = {} 
	pivotPairs = {}
	for e in reversed(edges):
		if(pivotcounter < 0 or pivots[pivotcounter].weight != e.weight or pivots[pivotcounter].polytop != e.polytop):
			faceList = getAllCofacets(e,dimension)
			columnV = []
			columnV.append(e)
			hq.heapify(faceList)
			while(True):
				while(faceList!=[]):
					pivot = faceList[0]
					hq.heappop(faceList)
					if(faceList!=[] and pivot.polytop==faceList[0].polytop):
						hq.heappop(faceList)
					else:
						hq.heappush(faceList,pivot)
						break
				if(faceList==[]):
					break
				elif pivot not in pivotPairs:
					pivotPairs[pivot] = e
					nextPivots.append(pivot)
					columnV.sort(key = letter_cmp_key)
					k =0
					for x in columnV:
						if(k+1 < len(columnV) and columnV[k]==columnV[k+1]):
							k+=1
						else:
							if e not in v:
								v[e] = [columnV[k]]
							else:
								f = v[e] + [columnV[k]]
								v[e] = f
							k+=1
					if(e.weight != pivot.weight):
						bettitableentry = [dimension,min(pivot.weight, e.weight),max(pivot.weight, e.weight)]
						bettieTable.append(bettitableentry)
					break
				else:
					if v and pivotPairs:
						for polytopnp in v[pivotPairs[pivot]]:
							columnV.append(polytopnp)
							faces =  getAllCofacets(polytopnp,dimension)
							for fc in faces:
								faceList.append(fc)
						hq.heapify(faceList)
		else:
			pivotcounter = pivotcounter - 1
	return nextPivots
	
def fastpersistance(polytopalcomplex):
	vertices = getDimEdges(0)
	edges = getDimEdges(1)
	pivots  = minimumspanningtree(edges)
	
	for d  in range(1,dim):
		
		if(d != 1):
			 edges = getDimEdges(d)
		pivots = persistenceByDimension(edges, pivots, d)
	return


def computedistancematrix(points):    #compute Distance Matrix
	minimumdistance = 99999
	maximumdistance = 0
	for x in range(0,len(points)):
		distance = []
		incidence = []
		for y in range(x,len(points)):
			di = math.dist(points[x],points[y])
			distance.append(di)
			incidence.append(0)
			if(di>maximumdistance):
				maximumdistance = di
			if(di<minimumdistance and di != 0):
				minimumdistance = di
		distancematrix.append(distance)
		delaunayincidence.append(incidence)
	return distancematrix

def distanceindex(x,y): # compute Distance by index b/w two points
	if(x>y):
		return distancematrix[y][x-y]
	else:
		return distancematrix[x][y-x]
def incidenceindexreturn(x,y): # compute Distance by index b/w two points
	if(x>y):
		return delaunayincidence[y][x-y]
	else:
		return delaunayincidence[x][y-x]
def setincidenceindex(x,y): # Update incedence by index b/w two points
	if(x>y):
		delaunayincidence[y][x-y] =1
	else:
		delaunayincidence[x][y-x] =1
		
def intersection(lst1, lst2):  #return intersection of two lists
    return list(set(lst1) & set(lst2))
 
def union(lst1,lst2):        #return union of two lists
	return list(set().union(lst1 , lst2))

def neighbours(polytop,simplices,dim1):    #report all the neighbouring polytopes that shares a face with polytope
	polytop = list(polytop)
	adjacency = []
	for y in simplices:
		y = list(y)
		if(polytop != y):
			if(len(intersection(polytop,y))==dim1):
				adjacency.append(y)
	return adjacency

def mergeneighbors(polytop,simplices,pointscoord):    #merge neigbours simplices that results in maximum convex polytop with neighbors
	neighbors = neighbours(polytop,simplices,len(pointscoord[0]))
	for x in neighbors:
		y = union(x,polytop)
		hull = ConvexHull([pointscoord[i] for i in y])
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


def optimalconvexization(simplices,pointscoord):   #Repeate the Maximum convexization for each simplex and returned the sorted list of convex polytopes by weight and vertices
	convexparts = []
	maxdist = []
	averagedist = []
	sizecd = []
	for i in range(0,len(simplices)):
		x= list(simplices[i])
		distance =0
		maxval = 0
		while True:
			x1 = list(mergeneighbors(x,simplices,pointscoord))
			if(x1==x):
				break
			x = x1
		if x not in convexparts:
			maxval=0
			distance = 0
			convexparts.append(x)
			hull = ConvexHull([pointscoord[i] for i in x])
			for p in hull.simplices:
				distance = distance + math.dist(pointscoord[x[p[0]]],pointscoord[x[p[1]]])
				if(maxval < math.dist(pointscoord[x[p[0]]],pointscoord[x[p[1]]])):
					maxval = math.dist(pointscoord[x[p[0]]],pointscoord[x[p[1]]])
			maxdist.append(maxval)
			averagedist.append(distance)
			sizecd.append(len(x))
		k=set()
		for tp in convexparts:
			for lm in tp:
				k.add(lm)
		if(len(k)==len(pointscoord)):
			break;
	df = pd.DataFrame(list(zip(convexparts,maxdist,averagedist,sizecd)),columns =['convexpart', 'maxdist','averagedist','sizecd'])
	df = df.sort_values(by = 'sizecd',ascending = False)
	df = df.sort_values(by = 'maxdist',ascending = True)
	return df

def iterativeconvexization(simplices,dim,pointscoord):   #Keep convex decomposition that minimizes convex polytopes required and minimizes maximum weight edge
	valid = []
	convexpartsunion = []
	simplices = list(simplices)
	remainingsimplices = simplices
	premaining = 0
	while(True):
		df = optimalconvexization(remainingsimplices,pointscoord)
		i=0
		pointsaddressed = []
		for x in df['convexpart'].tolist():
			check = 1
			for t in convexpartsunion:
				if(len(intersection(t,x)) > dim):
					check = 0
			if(check==1):
				valid.append(i)
				hull = ConvexHull([pointscoord[i] for i in x])
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
		if(remaining==premaining):
			if(remaining !=0):
				for r in remainingsimplices:
					convexpartsunion.append(r)
			break
		premaining = remaining
	convexpartssorted = []
	for x in convexpartsunion:
		convexpartssorted.append(sorted(x))
	return convexpartssorted


def fathestpointpolytopeindex(points,chebXc,pointscoord):
	farthestpoint = 0
	dist =0
	for i in range(0,len(points)):
		if(math.dist(pointscoord[i],chebXc)>dist):
			dist = math.dist(pointscoord[i],chebXc);
			farthestpoint = i
	return farthestpoint

def nearesttofarthestpoint(points,farthestpoint,pointscoord):
	nearesttofarthest = 0
	secdist =99999
	for i in range(0,len(points)):
		if(i!=farthestpoint):
			if(math.dist(pointscoord[i],pointscoord[farthestpoint])<secdist):
				secdist = math.dist(pointscoord[i],pointscoord[farthestpoint]);
				nearesttofarthest = i
	return nearesttofarthest	

def hyperplaneequation(chebXc,fpoint,ntofpoint): # Equation of the plane orthogonal to farthest points and passes through nearesttofarthest
	normal_vector = fpoint-chebXc
	constant=0
	for i in range(0,len(chebXc)):
		constant -= (normal_vector[i]*ntofpoint[i])
	return normal_vector,constant
	
def plt_sphere(c, r,ax):
    ax = fig.gca(projection='3d')

    # draw sphere
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
    x = r*np.cos(u)*np.sin(v)
    y = r*np.sin(u)*np.sin(v)
    z = r*np.cos(v)

    ax.plot_surface(x-c[0], y-c[1], z-c[2], color=np.random.choice(['g','b']), alpha=0.15*np.random.random())


def computesteriographicprojection(points,incenter,farthestpoint,nearesttofarthest,projectedcoordinates):
	pts = []
	cofficient,constant = hyperplaneequation(incenter,projectedcoordinates[farthestpoint],projectedcoordinates[nearesttofarthest])
	for i in range(0,len(points)):
		diff = projectedcoordinates[i]-projectedcoordinates[farthestpoint]
		constterm =0;
		diffterm = 0;
		for x in range(0,len(cofficient)):
			constterm+= projectedcoordinates[farthestpoint][x]*cofficient[x]
			diffterm+= diff[x]*cofficient[x]
		if(diffterm!=0):
			t = ((-1)*constant - constterm)/diffterm
			pnt = t*(projectedcoordinates[i]-projectedcoordinates[farthestpoint])+projectedcoordinates[farthestpoint]
			pts.append(pnt)
	return pts

def PCAprojection(pts):
	pca = PCA(n_components=len(pts[0])-1).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x
    # will need to add farthest point simplice

def generatesimplexfacets(poly,dimen):
	return list(combinations(poly, dimen))
	

def hyperplane(points):
   X=np.matrix(points)
   k=np.ones((len(points[0]),1))
   a=np.matrix.dot(np.linalg.inv(X), k)
   
   return a
   
def solutionlinalg(matrixA,matrixB):
	return np.linalg.inv(matrixA)*matrixB
	 

def findfarthestfacetprojection(polytop):
	polytopalpoints = [inputpoints[x] for x in polytop]
	
	polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopalpoints)
	polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
	chebR = polyhalfspace.chebR
	chebXc = polyhalfspace.chebXc
		
	faces = ConvexHull(polytopalpoints).simplices
	farthest = 0
	farthestface = []
	centroid1 = []
	for edges in faces:
		pts = [inputpoints[x] for x in edges]
		centroid = sum(np.transpose(list(list(x) for x in list(np.transpose(pts)))))/len(pts[0])
		if(math.dist(centroid,chebXc)>farthest):
			farthest = math.dist(centroid,chebXc)
			farthestface = edges
			centroid1 = centroid
	
	farthestface = list(farthestface)
	farthestfaceslist = []
	for x in farthestface:
		faceedge = [t for t in farthestface if t != x]
		for y in faces:
			if len(intersection(faceedge,y))==len(faceedge):
				if(not list(y) == list(farthestface)):
					farthestfaceslist.append(y)
	
	hyperplaneeq = []
	for x in farthestfaceslist:
         hyperplaneeq.append(hyperplane([inputpoints[t] for t in x]))
	matrixA = []
	matrixB = []
	for x in hyperplaneeq:
		matrixA.append([ g for t in np.array(x) for g in np.array(t)])
		matrixB.append([1])	
	matrixA = np.matrix(matrixA)
	matrixB = np.matrix(matrixB)
	projection_point = np.transpose(solutionlinalg(matrixA,matrixB))
	projection_point = [y for j in np.array(projection_point) for y in np.array(j)]
	
	stereoprojection_point = (centroid1+projection_point)/2
	stereoprojection_plane = hyperplane([inputpoints[t] for t in farthestface])
	cofficient = [y for j in np.array(stereoprojection_plane) for y in np.array(j)]
	constant  = 1
	'''
	p=0
	for i in range(0,len(polytop)):
		product = [a*b for a,b in zip(inputpoints[i],cofficient)]
		print(sum(product))
		if(sum(product)>.9):
			print(p)
			p=p+1
	'''
	pts = []
	for i in range(0,len(polytop)):
		diff = inputpoints[i]-stereoprojection_point
		constterm =0;
		diffterm = 0;
		for x in range(0,len(cofficient)):
			constterm+= stereoprojection_point[x]*cofficient[x]
			diffterm+= diff[x]*cofficient[x]
		if(diffterm!=0):
			t = ((-1)*constant - constterm)/diffterm
			pnt = t*(inputpoints[i]-stereoprojection_point)+stereoprojection_point
			pts.append(pnt)
	return PCAprojection(pts),stereoprojection_point
def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    points = np.transpose(points)
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return ctr, np.linalg.svd(M)[0][:,-1]

def projectonplane(point,planenormal,planecentroid):
	
	u = []
	for x,y in zip(point,planecentroid):
		u.append(x-y)
	n_norm = np.sqrt(sum(planenormal**2))    
	proj_of_u_on_n = (np.dot(u, planenormal)/n_norm**2)  
	return u - proj_of_u_on_n
	
def curvaturereducedconvexization(convexdecomposition,coordinates_points):
	finalizedpolytopes = []
	finalizedpolytopesspp = []
	convex_decomposition = convexdecomposition
	while len(convex_decomposition)!=0:
		convex_decomposition = sorted(convex_decomposition,key=len,reverse=True)
		y = convex_decomposition.pop(0)
		finalizedpolytopes.append(y)
		steriographic_projection ,projpoint = findfarthestfacetprojection(y)
		finalizedpolytopesspp.append(projpoint)
		faces = Delaunay(steriographic_projection).simplices
		cd = sorted(iterativeconvexization(faces,len(steriographic_projection[0]),steriographic_projection),key=len,reverse=True)
		th = 0
		for c in cd:
			ctr, normal  = planeFit([coordinates_points[k] for k in c])
			projectedpoints = []
			pts = []
			for t in c:
				projectedpoints.append(projectonplane(coordinates_points[t],normal,ctr))
				pts.append(t)
			th =0
			new_points = []
			faces1 = []
			for t in convex_decomposition:		
				if(t!=y):
					if(len(intersection(c,t))!=0):
						new_points = [r for r in t if r not in pts]
						for t1 in new_points:
							if t1 not in pts:
								pts.append(t1)
								projectedpoints.append(coordinates_points[t1])
						th = th +1
						faces1.append(t)
			facesindex = [x for x in range(0,len(pts))]
			new_faces = []
			for f in faces1:
				flist = []
				for fc in f:
					flist.append(pts.index(fc))
				new_faces.append(flist)
			cd1 = sorted(iterativeconvexization(new_faces,len(projectedpoints[0]),projectedpoints),key=len,reverse=True)
			mappedcd = []
			for f in cd1:
				flist = []
				for fc in f:
					flist.append(pts[fc])
				mappedcd.append(flist)
			updateconvexdecomposition = [tpr for tpr in convex_decomposition if tpr not in faces1]
			for tpr1 in mappedcd:
				updateconvexdecomposition.append(tpr1)
			convex_decomposition[:] = updateconvexdecomposition
	
	return finalizedpolytopes,finalizedpolytopesspp

'''
	
	reducedconvexization = []
	print(len(convexdecomposition))
	color = []
	n = 3
	for i in range(len(convexdecomposition)):
		color.append('#%06X' % randint(0, 0xFFFFFF))
	#tp = input()
	i=0
	fig = plt.figure(figsize=(10,7))                    #3
	ax=Axes3D(fig)                                      #3
	#ax.plot([i,i],[j,j],[k,h],color = 'g')
	#plt.show()
	X1 = np.transpose(inputpoints)[0]        #3
	Y1 = np.transpose(inputpoints)[1]        #3 
	Z1 = np.transpose(inputpoints)[2]        #3
	points = [x for x in range(0,len(inputpoints))]
	
	for x in convexdecomposition:
		#print(len(x))
		#print(x)
		if(len(x)>20):
			polytopal_points=np.array([inputpoints[y] for y in x])                                  #3
			weight = 0 #delauney retension
			points1 = [y for y in range(0,len(x))]
			polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopal_points)
			polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
			chebR = polyhalfspace.chebR
			chebXc = polyhalfspace.chebXc
			farthestpoint = fathestpointpolytopeindex(points1,chebXc,polytopal_points)
		
			nearesttofarthest = nearesttofarthestpoint(points1,farthestpoint,polytopal_points)
		
			steriographic_projection = computesteriographicprojection(points1 ,chebXc,farthestpoint,nearesttofarthest,polytopal_points)
			pca_x = PCAprojection(steriographic_projection)
	    
			hull = ConvexHull(polytopal_points)		
			for p in hull.simplices:
				if(len(p)>3):
					print(p)
				for t in range(0,len(p)):
					for e in range(t+1,len(p)): 
						#plt.plot([points[x[p[0]]][0],points[x[p[1]]][0]],[points[x[p[0]]][1],points[x[p[1]]][1]])
						ax.plot([X1[x[p[t]]],X1[x[p[e]]]], [Y1[x[p[t]]],Y1[x[p[e]]]], [Z1[x[p[t]]],Z1[x[p[e]]]],alpha=0.7,color=color[i])
			i=i+1
			print(i," ",len(x))		
			#plt.scatter(points[:,0],points[:,1],s=10,color='black')
			ax.scatter(np.transpose(polytopal_points)[0], np.transpose(polytopal_points)[1], np.transpose(polytopal_points)[2])
#			plt.show()
			it= input()
			faces1 = []
			boundaryfaces = []
			for R in hull.simplices:
				if farthestpoint not in R:
					faces1.append(R)
				else:
					boundaryfaces.append(R)
			faces = []
			for c in faces1:
				fc = []
				for k in c:
					if(k<farthestpoint):
						fc.append(k)
					else:
						fc.append(k-1)
				faces.append(fc)
			convexdecomposedfaces12 = iterativeconvexization(faces,len(pca_x[0]),pca_x)
			for p in convexdecomposedfaces12:
				if(len(p)>3):
					print(p)
				pts = [inputpoints[x[t]] for t in p]
				#for t in range(0,len(p)):
					#for e in range(t+1,len(p)): 
					#	plt.plot([points[x[p[0]]][0],points[x[p[1]]][0]],[points[x[p[0]]][1],points[x[p[1]]][1]])
					#	#ax.plot([X1[x[p[t]]],X1[x[p[e]]]], [Y1[x[p[t]]],Y1[x[p[e]]]], [Z1[x[p[t]]],Z1[x[p[e]]]],alpha=0.1,color=color[i])
			#plt.scatter(points[:,0],points[:,1],s=10,color='black')
			#ax.scatter(np.transpose(polytopal_points)[0], np.transpose(polytopal_points)[1], np.transpose(polytopal_points)[2])
				ax.scatter(np.transpose(pts)[0], np.transpose(pts)[1], np.transpose(pts)[2])
			plt.show()
			it = input()
'''            
class DisjointSet:
    def __init__(self, n):
        self.parent = [i for i in range(n)]
        self.rank = [1 for _ in range(n)]
    
    # make a and b part of the same component
    # union by rank optimization
    def union(self, a, b):
        pa = self.find(a)
        pb = self.find(b)
        if pa == pb: return
        if self.rank[pa] > self.rank[pb]:
            self.parent[pb] = pa
            self.rank[pa] += self.rank[pb]
        else:
            self.parent[pa] = pb
            self.rank[pb] += self.rank[pa]
    
    # find the representative of the 
    # path compression optimization
    def find(self, a):
        if self.parent[a] == a:
            return a
        
        self.parent[a] = self.find(self.parent[a])
        return self.parent[a]

class orderedarraylistnode(object):
	def __init__(self,polytop,weight):
		self.polytop = sorted(polytop,reverse=False)
		self.weight = weight
	
	def __hash__(self):
		return hash((tuple(self.polytop), self.weight))
	
	def __getitem__(self, item):
		return [self.polytop,self.weight]
	
	def __lt__(self, nxt):
		if self.weight < nxt.weight:
			return True
		elif self.weight == nxt.weight:
			for (x,y) in zip(reversed(self.polytop),reversed(nxt.polytop)):
				if(x > y):
					return True
		return False
	def __eq__(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		return self.polytop == other.polytop and self.weight == other.weight
        
        
class polytope:
	def __init__(self,polytope1,polyunmapped,weight1,farthestpoint1,convexdecomposedfaces,mappedconvexdecomposedfaces,projpoints,dimension,pca_x):
		self.polytopevertices = polytope1  # list of coordinates of polytope vertices
		self.polytopunmappedvertices = polyunmapped
		self.weight = weight1  #weight of polytope
		self.farthestpoint = farthestpoint1 #polytope fathest point from incenter coordinates
		self.unmappedconvexdecomposed = convexdecomposedfaces
		self.convexdecomposedfaces =  mappedconvexdecomposedfaces  #PCA Dimensionality Reduction for lower order faces
		self.projectedpoints = projpoints
		self.dimension = dimension # vertices in highest simplex
		self.childrens = []
		self.pca_x = pca_x

class tree:
	def createheadpolytope(self,points):
		poly = points
		weight = 0
		farthestpoint = 9999999
		faces =  Delaunay(inputpoints).simplices
		for x in faces:
			for k in generatesimplexfacets(x,2):
				setincidenceindex(k[0],k[1])
		convexdecomposedfaces1 = iterativeconvexization(faces,len(inputpoints[0]),inputpoints)
		convexdecomposedfaces,projpoints = curvaturereducedconvexization(convexdecomposedfaces1,inputpoints)
		print(len(convexdecomposedfaces))
		print(len(convexdecomposedfaces1))
		print(len(projpoints))
		i=0
		for x in sorted(convexdecomposedfaces,key=len,reverse=True):
			print(x)
			print(i)
			i=i+1
		t = input()
		dimension = len(inputpoints[0])
		pca_x = inputpoints
		return polytope(poly,poly,weight,farthestpoint,convexdecomposedfaces,convexdecomposedfaces,projpoints,dimension,pca_x)

	def createpolytope(self,points,points1,pca_x1,farthest,dimension):
		poly = points
		polyunmapped = points1

		polytopal_points = []
		if(farthest == -1 or len(points1) == dimension+1):
			mappedconvexfaces = generatesimplexfacets(poly, dimension)
			convexfaces = generatesimplexfacets([i for i in range(0,len(points))], dimension)
			return polytope(poly,points1,0,-1,convexfaces,mappedconvexfaces,dimension-1,[])
			
		for i in points1:
			polytopal_points.append(pca_x1[i])
		weight = 0 #delauney retension
		polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopal_points)
		polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
		chebR = polyhalfspace.chebR
		chebXc = polyhalfspace.chebXc
		farthestpoint = fathestpointpolytopeindex(points1,chebXc,polytopal_points)
		
		nearesttofarthest = nearesttofarthestpoint(points1,farthestpoint,polytopal_points)
		
		steriographic_projection = computesteriographicprojection(points1 ,chebXc,farthestpoint,nearesttofarthest,polytopal_points)
		pca_x = PCAprojection(steriographic_projection)
	    
		hull = ConvexHull(polytopal_points).simplices
		faces1 = []
		boundaryfaces = []
		for x in hull:
			if farthestpoint not in x:
				faces1.append(x)
			else:
				boundaryfaces.append(x)
		faces = []
		for c in faces1:
			fc = []
			for k in c:
				if(k<farthestpoint):
					fc.append(k)
				else:
					fc.append(k-1)
			faces.append(fc)
		if(len(pca_x[0]) <= 1):
			convexdecomposedfaces = faces		
		else:
			convexdecomposedfaces1 = iterativeconvexization(faces,len(pca_x[0]),pca_x)
			convexdecomposedfaces = curvaturereducedconvexization(convexdecomposedfaces1,pca_x)
		mappedconvexedfaces = []
		for a in convexdecomposedfaces:
			list1 = []
			for y in a:
				if(y>=farthestpoint):
					list1.append(points[y+1])
				else:
					list1.append(points[y])
			mappedconvexedfaces.append(sorted(list1))
		for x in boundaryfaces:
			list1 = []
			list2 = []
			for j in x:
				list1.append(points[j])
				list2.append(j)	
			convexdecomposedfaces.append(sorted(list2))
			mappedconvexedfaces.append(sorted(list1))

		dimension = len(pca_x[0])
		return polytope(poly,polyunmapped,weight,farthestpoint,convexdecomposedfaces,mappedconvexedfaces,dimension,pca_x)
	
	def createbootomup(self,poly):  #needs to be implemented
		i = 0
		queue = []
		queue.append(poly)
		polytopaltree = []
		polytopallevel = []
		xdimension = []	  
		currentdim = dim
		queuelevel = []
		level = 0
		previouslevel = 0
		queuelevel.append(level)
		k=0
		while(queue !=[] and  currentdim>0):
			poly = queue.pop(0)
			level1 = queuelevel.pop(0)
			level = level1
			currentdim = poly.dimension
			k = 0
			for polytop in poly.convexdecomposedfaces:
				if(len(polytop) >=1):
					if(len(polytop)==poly.dimension+1 or	 len(polytop) ==1):
						mappedconvexfaces = generatesimplexfacets(polytop,poly.dimension)
						convexfaces = generatesimplexfacets([i for i in range(0,len(polytop))],poly.dimension)
						newchild = polytope(polytop,poly.unmappedconvexdecomposed[k],0,-1,convexfaces,mappedconvexfaces,poly.dimension-1,[])
					else:
						newchild = self.createpolytope(polytop,poly.unmappedconvexdecomposed[k],poly.pca_x,poly.farthestpoint,poly.dimension-1)
					i=i+1
					if(i%10000==0):
						print(i)
					poly.childrens.append(newchild)
					queue.append(newchild)
					queuelevel.append(level+1)
					tuple1 = ()
					if(level>previouslevel):
						res = []
						res1 = []
						for x in xdimension:
							if(type(x)!=type(xdimension)):
								if(type(x) == type(tuple1)):
									res1.append(list(x))
								else:
									res1.append(x.tolist())
							else:
								res1.append(x)
						[res.append(x) for x in res1 if x not in res]
						polytopaltree.append(res)
						xdimension = []
						previouslevel = level
					xdimension.append(polytop)
					k=k+1
		vertices = []
		for i in range(0,datasize):
			vertices.append([i])
		polytopaltree.append(vertices)    
		return polytopaltree
	def display(self,root):  # needs to be implimented
		print("polytope vertices")
		print(root.polytopevertices)
		print("polytope Weight")
		print(root.weight)
		print("Half Spaces")
		print(root.polytophalyspaces)
		print("Polytope Equation")
		print(root.polyhalfspace)
		print("polytope inRadius")
		print(root.chebR)
		print("polytop incenter")
		print(root.chebXc)
		print("polytope farthest point")
		print(root.farthestpoint)
		print("polytope Stereographic projections")
		print(root.steriographic_projection)
		print("polytope faces")
		print(root.faces)
		print("convex decomposed faces")
		print(root.convexdecomposedfaces)
		print("Dimensions")
		print(root.dimension)
		return

distancematrix = []	  #initializae distance matrix
delaunayincidence = []
head = tree()  # initialize a polytopal tree
#read file
#print("Here")
'''
inputpoints = []
with open("lion200.csv", "r") as f:
	reader = csv.reader(f)
	for line in reader:
		k=[]
		for x in line:
			k.append(float(x))
		inputpoints.append(k)
			
#print(inputpoints)

inputpoints = []
def genPermutahedron(a, size):
    # if size becomes 1 then prints the obtained
    # permutation
    if size == 1:
        inputpoints.append(np.array(a))
        return
 
    for i in range(size):
        genPermutahedron(a, size-1)
 
        # if size is odd, swap 0th i.e (first)
        # and (size-1)th i.e (last) element
        # else If size is even, swap ith
        # and (size-1)th i.e (last) element
        if size & 1:
            a[0], a[size-1] = a[size-1], a[0]
        else:
            a[i], a[size-1] = a[size-1], a[i]
            
import itertools

def generatehypercube(dimensions):
    return list(itertools.product('10', repeat=dimensions))

b = [0,1,2,4]
size = len(b)
genPermutahedron(b,size)
#dim =4
#b=([1,1,1,1])
#inputpoints = [[int(y) for y in x] for x in generatehypercube(dim)]
#print(inputpoints)
#ty = input()
x, y = np.meshgrid(np.linspace(-1,1,len(b)-1), np.linspace(-1,1,24))
dst = np.sqrt(x*x+y*y)

sigma = 1
muu = 0.000
 
# Calculating Gaussian array
gauss = np.exp(-( (dst-muu)**2 / ( 2.0 * sigma**2 ) ) )

pca_x = PCAprojection(inputpoints)
inputpoints = pca_x
#print(inputpoints)
#print(gauss)

inputpoints = np.array(inputpoints) + gauss
'''
#print(inputpoints)
#print(inputpoints)            
#dim= len(inputpoints[0])  #dimension of a data
dim = 4
#datasize = len(inputpoints)
datasize = 40
inputpoints = tadasets.dsphere(n=datasize, d=dim-1, r=1, noise=1)  # generate d-sphere datatset
#inputpoints = np.array([[3,2,7],[5,4,3],[-3,-1,6],[-2,-3,4],[-3,2,6],[4,-4,3],[2,-2,5],[0,0,-8],[-5,5,1],[-5,-5,1]])
#inputpoints = np.array([[-0.33503,-0.93161,0.14092],[-0.58394,0.6616,0.47041],[-0.60231,-0.10734,-0.79101],[0.5778,0.58913,0.56486],[0.96496,0.24644,-0.09013],[-0.49422,-0.0695,0.86655],[0.9337,0.34562,-0.09359],[0.95858,-0.19768,-0.20505],[-0.75412,-0.62654,-0.19683],[-0.38701,0.41265,0.82459],[0.86506,0.40221,0.29983],[0.69299,0.37824,0.61376],[-0.35164,-0.6592,-0.66468],[-0.42816,0.16204,-0.88906],[0.93148,-0.17783, 0.31738],[ 0.06839, 0.89008 ,0.45064],[ 0.3142,  0.7662 ,-0.56055],[ 0.89106,0.07213,0.44813],[-0.43599,0.70498,0.55939],[-0.86108,0.49248,0.1265]])
points = [i for i in range(0,len(inputpoints))]
computedistancematrix(inputpoints)	#Populate distance matrix
dimension1 =dim
polytopal = True
simplicial = True
if(polytopal==True):
	root = head.createheadpolytope(points)  # if data lies in n-dimensions; this will make a surface of a lifted (n+1)-dimensional parbolied its projection will be delaunay triangulation
	polytopalarraylist  = head.createbootomup(root) 
	weights = []
	for dimensionlist in polytopalarraylist:
		distancedime = []
		for polytope in dimensionlist:
			maxedge = 0
			if(len(polytope)!= 1):
				edges = generatesimplexfacets(polytope,2)
				for edge in edges:
					if(incidenceindexreturn(edge[0],edge[1])==1):
						if(maxedge<distanceindex(edge[0],edge[1])):
							maxedge = distanceindex(edge[0],edge[1])
				if(maxedge==0):
					distancedime.append(distanceindex(edge[0],edge[1]))#distanceindex(edge[0],edge[1]))
				else:
					distancedime.append(maxedge)
			else:
				distancedime.append(0)
		weights.append(distancedime)

        
        
	orderedpolytoparraylist = []        
	for (polytop, weight) in zip(polytopalarraylist, weights):
		dimensionwiseorderedpolytopes = set()
		for (poly,weig) in zip(polytop,weight):
			node = orderedarraylistnode(poly,weig)
			dimensionwiseorderedpolytopes.add(node)
		sorted_list = list(dimensionwiseorderedpolytopes)
		sorted_list.sort(key = letter_cmp_key)
		orderedpolytoparraylist.append(sorted_list)
		dimension1 = dim
	for x in orderedpolytoparraylist:
		print("polytopes of Dimension ::",dimension1,"  ::  ",len(x))
		#for p in x:
		#	print(p.polytop,p.weight)
		dimension1 = dimension1-1
	
	bettieTable = []
	fastpersistance(orderedpolytoparraylist)
	for x in bettieTable:
		print(x)
if(simplicial == True):
	orderedpolytoparraylist1 = []        
	dimension1 = 3

	for d in range(0,dimension1+1):
		dimensionwiseorderedpolytopes = set()
		if(d==0):
			for i in range(0,datasize):
				node = orderedarraylistnode([i],0)
				dimensionwiseorderedpolytopes.add(node)
		else:
			for i in orderedpolytoparraylist1[d-1]:
				for j in range(max(i.polytop)+1,datasize):
					weight1 = i.weight
					for p in i.polytop:
						if(distanceindex(p,j)>weight1):
							weight1 =  distanceindex(p,j)
					node = orderedarraylistnode(i.polytop+[j],weight1)
					dimensionwiseorderedpolytopes.add(node)
		sorted_list = list(dimensionwiseorderedpolytopes)
		sorted_list.sort(key = letter_cmp_key)
		orderedpolytoparraylist1.append(sorted_list)

	orderedpolytoparraylist = []        
	for x in range(dim,-1,-1):
		orderedpolytoparraylist.append(orderedpolytoparraylist1[x])
	
	#simplicial Complex
	for x in orderedpolytoparraylist:
		print("simplicies of Dimension ::",dimension1,"  ::  ",len(x))
		#for p in x:
		#	print(p.polytop,p.weight)
		dimension1 = dimension1-1

	bettieTable = []
	fastpersistance(orderedpolytoparraylist)
	for x in bettieTable:
		print(x)
