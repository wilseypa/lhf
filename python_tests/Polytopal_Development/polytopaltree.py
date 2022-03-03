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


class polytope:
	def __init__(self,polytope1,polyunmapped,weight1,farthestpoint1,convexdecomposedfaces,mappedconvexdecomposedfaces,dimension,pca_x):
		self.polytopevertices = polytope1  # list of coordinates of polytope vertices
		self.polytopunmappedvertices = polyunmapped
		self.weight = weight1  #weight of polytope
		self.farthestpoint = farthestpoint1 #polytope fathest point from incenter coordinates
		self.convexdecomposed = convexdecomposedfaces
		self.convexdecomposedfaces =  mappedconvexdecomposedfaces  #PCA Dimensionality Reduction for lower order faces
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
		convexdecomposedfaces = iterativeconvexization(faces,len(inputpoints[0]),inputpoints)
		dimension = len(inputpoints[0])
		pca_x = inputpoints
		return polytope(poly,poly,weight,farthestpoint,convexdecomposedfaces,convexdecomposedfaces,dimension,pca_x)

	def createpolytope(self,points,points1,pca_x1,farthest,dimension):
		poly = points
		polyunmapped = points1
		polytopal_points = []
		if(farthest ==-1 or len(points1) == dimension+1):
			mappedconvexfaces = generatesimplexfacets(poly, dimension)
			convexfaces = generatesimplexfacets([i for i in range(0,len(points))], dimension)
			return polytope(poly,points1,0,-1,convexfaces,mappedconvexfaces,dimension-1,[])
			
		for i in points1:
			if(len(pca_x1[0]) == dim):
				polytopal_points.append(pca_x1[i])
			elif(i>=farthest):
				polytopal_points.append(pca_x1[i-1])
			else:
				polytopal_points.append(pca_x1[i])
		weight = 0 #delauney retension
		polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopal_points)
		polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
		chebR = polyhalfspace.chebR
		chebXc = polyhalfspace.chebXc
		farthestpoint = fathestpointpolytopeindex(points1,chebXc,polytopal_points)
		#print(points[farthestpoint])
		
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
				if(k>=farthestpoint):
					fc.append(k-1)
				else:
					fc.append(k)
			faces.append(fc)
		if(len(pca_x[0]) <= 1):
			convexdecomposedfaces = faces		
		else:
			convexdecomposedfaces = iterativeconvexization(faces,len(pca_x[0]),pca_x)
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
					if(len(polytop)==poly.dimension or	 len(polytop) ==1):
						mappedconvexfaces = generatesimplexfacets(polytop,poly.dimension)
						convexfaces = generatesimplexfacets([i for i in range(0,len(polytop))],poly.dimension)
						newchild = polytope(polytop,poly.convexdecomposed[k],0,-1,convexfaces,mappedconvexfaces,poly.dimension-1,[])
					else:
						newchild = self.createpolytope(polytop,poly.convexdecomposed[k],poly.pca_x,poly.farthestpoint,poly.dimension)
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
dim= 3  #dimension of a data
datasize = 10
inputpoints = tadasets.dsphere(n=datasize, d=dim-1, r=1, noise=0)  # generate d-sphere datatset
inputpoints = np.array([[3,2,7],[5,4,3],[-3,-1,6],[-2,-3,4],[-3,2,6],[4,-4,3],[2,-2,5],[0,0,-8],[-5,5,1],[-5,-5,1]])
points = [i for i in range(0,len(inputpoints))]
computedistancematrix(inputpoints)	#Populate distance matrix
root = head.createheadpolytope(points)  # if data lies in n-dimensions; this will make a surface of a lifted (n+1)-dimensional parbolied its projection will be delaunay triangulation
polytopalarraylist  = head.createbootomup(root)  #will need to impliment this function to gro polytopal tree bottom up
dimension1 =dim
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
				distancedime.append(maxedge)#distanceindex(edge[0],edge[1]))
			else:
				distancedime.append(maxedge)
		else:
			distancedime.append(0)
	weights.append(distancedime)

class orderedarraylistnode(object):
	def __init__(self,polytop,weight):
		self.polytop = sorted(polytop,reverse=True)
		self.weight = weight
	
	def __hash__(self):
		return hash((tuple(self.polytop), self.weight))
	
	def __getitem__(self, item):
		return [self.polytop,self.weight]
	
	def __lt__(self, nxt):
		if self.weight > nxt.weight:
			return True
		elif self.polytop > nxt.polytop:
			return True
		else:
			return False
	
	def __eq__(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		return self.polytop == other.polytop and self.weight == other.weight
        
        
        
orderedpolytoparraylist = []        


for (polytop, weight) in zip(polytopalarraylist, weights):
	dimensionwiseorderedpolytopes = set()
	for (poly,weig) in zip(polytop,weight):
		node = orderedarraylistnode(poly,weig)
		dimensionwiseorderedpolytopes.add(node)
	sorted_list = sorted(dimensionwiseorderedpolytopes, key = lambda x: (x.weight, x.polytop))
	orderedpolytoparraylist.append(sorted_list)
	
for x in orderedpolytoparraylist:
	print("polytopes of Dimension ::",dimension1,"  ::  ")
	for p in x:
		print(p.polytop,p.weight)
	dimension1 = dimension1-1

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
    
bettieTable = []
def minimumspanningtree(edges1):
	# using Kruskal's algorithm to find the cost of Minimum Spanning Tree
	res = 0
	pivots = []
	mstSize = 0
	ds = DisjointSet(datasize)
	for x in edges1:
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

def getAllCofacets(polytop,dimension):
	returnlist = []
	for x in orderedpolytoparraylist[dim-dimension-1]:
		if(len(intersection(x.polytop,polytop)) == len(polytop)):
			returnlist.append(x)
	return returnlist
				
	
def persistenceByDimension( edges, pivots, dimension):
	pivotcounter = 0
	nextPivots = []	
	v = {} 
	#for x in pivots:
	#	print(x.polytop,x.weight)
	pivotPairs = {}	
	for e in edges:
		if(pivotcounter >= len(pivots) or pivots[pivotcounter].weight != e.weight or pivots[pivotcounter].polytop != e.polytop):
			faceList = getAllCofacets(e.polytop,dimension)
			#print(e.polytop)
			#for x in faceList:
			#	print(x.polytop,x.weight)
			#x=input()
			columnV = []
			columnV.append(e)
			hq.heapify(faceList)
			while(True):
				while(faceList!=[]):
					pivot = faceList[0]
					hq.heappop(faceList)
					if(faceList!=[] and pivot==faceList[0]):
						hq.heappop(faceList)
					else:
						hq.heappush(faceList,pivot)
						break
				if(faceList==[]):
					break
				elif pivot not in pivotPairs:
					pivotPairs[pivot] = e
					nextPivots.append(pivot)
					for x in columnV:
						print(x.polytop,x.weight)
					
					columnV = sorted(columnV, key = lambda x: (x.weight, x.polytop))
					for x in columnV:
						print(x.polytop,x.weight)
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
						print(len(v))
						print(pivotPairs[pivot])
						for polytopnp in v[pivotPairs[pivot]]:
							columnV.append(polytopnp)
							faces =  getAllCofacets(polytopnp.polytop,dimension)
							for fc in faces:
								print("k:",k)
								k =k+1
								faceList.append(fc)
						hq.heapify(faceList)
			print("here")
		else:
			pivotcounter+=1
	return nextPivots
	
def fastpersistance(polytopalcomplex):
	vertices = getDimEdges(0)
	edges = getDimEdges(1)
	pivots  = minimumspanningtree(edges)
	
	for d  in range(1,dim):
		if(d != 1):
			 edges = getDimEdges(d)

		pivots = persistenceByDimension(edges, pivots, d)
		for x in pivots:
			print(x.polytop,x.weight)
	return

fastpersistance(orderedpolytoparraylist)

for x in bettieTable:
	print(x)

