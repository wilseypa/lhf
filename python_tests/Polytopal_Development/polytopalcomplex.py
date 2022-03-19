from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import tadasets
import math
import numpy as np
import copy
from itertools import combinations
import pandas as pd
import pypoman.duality as polyfun
import polytope as pc
from sklearn.decomposition import PCA
from functools import cmp_to_key
import heapq as hq


distancematrix = []	  #initializae distance matrix
delaunayincidence = []

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


def generatesimplexfacets(poly,dimen):  #generate simplices of dimension dimen
	tp = list(combinations(poly, dimen))
	listreturn = []
	for x in tp:
		listreturn.append(list(x))
	return listreturn

def mapindices(fromindices,toindices):
	newindices = []
	for x in fromindices:
		lst = []
		for y in x:
			lst.append(toindices[y])
		newindices.append(lst)
	return np.array(newindices,dtype=object)
	
def hyperplane(points):
   X=np.matrix(points)
   k=np.ones((len(points[0]),1))
   a=np.matrix.dot(np.linalg.inv(X), k)
   
   return a
 
def solutionlinalg(matrixA,matrixB):
	return np.linalg.inv(matrixA)*matrixB  

def PCAprojection(pts):
	pca = PCA(n_components=len(pts[0])-1).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x
def findfarthestfacetprojection(coordinates_points):  #compute farthest facet steriographic projection
	polytopalpoints = coordinates_points
	
	polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopalpoints)
	polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
	chebR = polyhalfspace.chebR
	chebXc = polyhalfspace.chebXc
		
	faces = ConvexHull(polytopalpoints).simplices
	farthest = 0
	farthestface = []
	centroid1 = []
	for edges in faces:
		pts = [coordinates_points[x] for x in edges]
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
         hyperplaneeq.append(hyperplane([coordinates_points[t] for t in x]))
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
	stereoprojection_plane = hyperplane([coordinates_points[t] for t in farthestface])
	cofficient = [y for j in np.array(stereoprojection_plane) for y in np.array(j)]
	constant  = 1

	pts = []
	for i in range(0,len(coordinates_points)):
		diff = coordinates_points[i]-stereoprojection_point
		constterm =0;
		diffterm = 0;
		for x in range(0,len(cofficient)):
			constterm+= stereoprojection_point[x]*cofficient[x]
			diffterm+= diff[x]*cofficient[x]
		if(diffterm!=0):
			t = ((-1)*constant - constterm)/diffterm
			pnt = t*(coordinates_points[i]-stereoprojection_point)+stereoprojection_point
			pts.append(pnt)
	return PCAprojection(pts)

def generatelowerorderdecompositions(indices,coordinates):
	level = []
	indices = indices
	coords = coordinates
	points = []
	[points.append(x) for z in indices for x in z if x not in points]
	points = sorted(points)
	triangulation = ConvexHull(coords).simplices
	projectedpts = findfarthestfacetprojection(coords)
	convexdecomposition = iterativeconvexization(triangulation,len(projectedpts[0]),projectedpts)
	convexdecompositionmapped = mapindices(convexdecomposition,points)
	for (x,y) in zip(convexdecomposition,convexdecompositionmapped):
		x = sorted(x)
		lowercoordinates = [projectedpts[i] for i in x]
		decomposition = Delaunay(lowercoordinates).simplices
		mappedindices = mapindices(decomposition,y)
		level.append([mappedindices.tolist(),lowercoordinates])
	return level

def addEdgesAndVertices(polytopaltree):
	level = []
	for x in polytopaltree[-1]:
		indices = x[0]
		coords = x[1]
		points = []
		[points.append(x) for z in indices for x in z if x not in points]
		points = sorted(points)
		if(len(points)>3):
			triangulation = ConvexHull(coords).simplices
			mappedindices = mapindices(triangulation,points)
			[level.append(sorted(t)) for t in mappedindices.tolist() if sorted(t) not in level]
		else:
			facets = generatesimplexfacets(points,2)
			[level.append(sorted(t)) for t in facets if sorted(t) not in level]
	polytopaltree.append(level)
	print("level::",1)
	level = []
	for x in range(0,len(inputpoints)):
		level.append([x])
	polytopaltree.append(level)
	print("level::",0)
	return polytopaltree

def createnormalpolytopaltree(inputpoints):
	computedistancematrix(inputpoints)
	coordinates = []
	coordinates = copy.deepcopy(inputpoints)
	polytopaltree = []
	#initialize the top level polytopes
	triangulation =  Delaunay(coordinates).simplices
	for x in triangulation:
		for k in generatesimplexfacets(x,2):
			setincidenceindex(k[0],k[1])
	convexdecomposition = iterativeconvexization(triangulation,len(coordinates[0]),coordinates)
	level = []
	for x in convexdecomposition:
		x = sorted(x)
		lowercoordinates = [coordinates[i] for i in x]
		decomposition = Delaunay(lowercoordinates).simplices
		mappedindices = mapindices(decomposition,x)
		level.append([mappedindices.tolist(),lowercoordinates])
	polytopaltree.append(level)
	print("level::",dim)
	#cascade down to Triangles
	for d in range(1,dim-1):
		convexdecompositions = copy.deepcopy(polytopaltree[d-1])
		level = []
		levelcheck = []
		lenthr = 0
		for cd in convexdecompositions:
			print(cd)
			if(len(cd[0])!=1):
				lowercd = ConvexHull(cd[1]).simplices
				print(lowercd)
				for x in lowercd:
					y = sorted(x)
					if y not in levelcheck:
						levelcheck.append(y)
						level.append([[y],[]])
				lenthr = lenthr + len(lowercd)
			else:
				facets = generatesimplexfacets(cd[0][0],len(cd[0][0])-1)
				for x in facets:
					y = sorted(x)
					y= [y]
					if y not in levelcheck:
						levelcheck.append(y)
						level.append([y,[]])
				lenthr = lenthr + len(facets)
		polytopaltree.append(level)
		print("level ::",(dim-d))
	#Add edges and vertices
	polytopaltree = addEdgesAndVertices(polytopaltree)
	return polytopaltree
def createpolytopaltree(inputpoints):
	computedistancematrix(inputpoints)
	coordinates = []
	coordinates = copy.deepcopy(inputpoints)
	polytopaltree = []
	#initialize the top level polytopes
	triangulation =  Delaunay(coordinates).simplices
	for x in triangulation:
		for k in generatesimplexfacets(x,2):
			setincidenceindex(k[0],k[1])
	convexdecomposition = iterativeconvexization(triangulation,len(coordinates[0]),coordinates)
	level = []
	for x in convexdecomposition:
		x = sorted(x)
		lowercoordinates = [coordinates[i] for i in x]
		decomposition = Delaunay(lowercoordinates).simplices
		mappedindices = mapindices(decomposition,x)
		level.append([mappedindices.tolist(),lowercoordinates])
	polytopaltree.append(level)
	print("level::",dim)
	#cascade down to Triangles
	for d in range(1,dim-1):
		convexdecompositions = copy.deepcopy(polytopaltree[d-1])
		level = []
		levelcheck = []
		lenthr = 0
		for cd in convexdecompositions:
			if(len(cd[0])!=1):
				lowercd = generatelowerorderdecompositions(cd[0],cd[1])
				for x in lowercd:
					y = sorted(x[0])
					if y not in levelcheck:
						levelcheck.append(y)
						level.append(x)
				lenthr = lenthr + len(lowercd)
			else:
				facets = generatesimplexfacets(cd[0][0],len(cd[0][0])-1)
				for x in facets:
					y = sorted(x)
					y= [y]
					if y not in levelcheck:
						levelcheck.append(y)
						level.append([y,[]])
				lenthr = lenthr + len(facets)
		polytopaltree.append(level)
		print("level ::",(dim-d))
	#Add edges and vertices
	polytopaltree = addEdgesAndVertices(polytopaltree)
	return polytopaltree
	
class orderedarraylistnode(object):
	def __init__(self,polytop,weight):
		self.polytopparts = polytop
		polytop1 = []
		if(type(polytop[0][0])!=type(np.int32(1)) and type(polytop[0][0])!=type(1)):
			[polytop1.append(x) for y in polytop for x in y[0] if x not in polytop1]
		else:
			[polytop1.append(x) for x in polytop[0] if x not in polytop1]
		polytop1 = sorted(polytop1)
		self.polytop = sorted(polytop1,reverse=False)
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

class orderedarraylistnodeorig(object):
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
        
def assignweights(polytopaltree):
	orderedpolytoparraylist = []
	for dimensionlist in polytopaltree:
		dimensionwiseorderedpolytopes = set()
		for polytope in dimensionlist:
			if(type(polytope[0])!=type(1)):
				if(type(polytope[0])!=type(np.int32(1))):
					polytp = []
					maxedge = 0
					for poly in polytope[0]:
						localmax =0
						edges = generatesimplexfacets(poly,2)
						for edge in edges:
							if(incidenceindexreturn(edge[0],edge[1])==1):
								if(localmax<distanceindex(edge[0],edge[1])):
									localmax = distanceindex(edge[0],edge[1])
						polytp.append([poly,localmax])
						if(maxedge<localmax):
							maxedge = localmax
					node = orderedarraylistnode(polytp,maxedge)
				else:
					maxedge = distanceindex(polytope[0],polytope[1])
					node = orderedarraylistnode([polytope,maxedge],maxedge)
			else:
				node = orderedarraylistnode([polytope,0],0)
			dimensionwiseorderedpolytopes.add(node)
		sorted_list = list(dimensionwiseorderedpolytopes)
		sorted_list.sort(key = letter_cmp_key)
		orderedpolytoparraylist.append(sorted_list)
	return orderedpolytoparraylist


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
			bettieTable.add(tuple(bettitableentry))
			if(mstSize >= len(inputpoints)-1):
				for i in range(0,len(inputpoints)):
					if(ds.find(i) == i):
						bettitableentry = [0,0,'inf']
						bettieTable.add(tuple(bettitableentry))
				return pivots
	return pivots

def getDimEdges(dimension):
	edg =  polytopaltree[dim -dimension]
	return edg
def getAllCofacets(e,dimension):    #will update this function a little
	returnlist = []
	for x in polytopaltree[dim-dimension-1]:
		if(len(intersection(x.polytop,e.polytop)) == len(e.polytop)):
			#if(len(x.polytop) == dimension+2):
			returnlist.append(x)
	return returnlist

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
			i = 0
			while(True):
				print(i)
				i=i+1
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
					if(e.weight < pivot.weight):
						bettitableentry = [dimension,min(pivot.weight, e.weight),max(pivot.weight, e.weight)]
						bettieTable.add(tuple(bettitableentry))
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
def updatedpivots(pivots):
	print(pivots)
	print(type(pivots[0]))
	k = input()
def fastpersistance(polytopalcomplex):
	vertices = getDimEdges(0)
	edges = getDimEdges(1)
	pivots  = minimumspanningtree(edges)
	
	for d  in range(1,dim):
		
		if(d != 1):
			 edges = getDimEdges(d)
		pivots = persistenceByDimension(edges, pivots, d)
		#pivots = updatedpivots(pivots)
	return
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
def disintegrate(polytopaltree):
	orderedpolytoparraylist = []
	for dimensionlist in polytopaltree:
		dimensionwiseorderedpolytopes = set()
		for polytope in dimensionlist:
			if(type(polytope.polytopparts[0][0]) == type([])):
				for part in polytope.polytopparts:
					node = orderedarraylistnodeorig(part[0],part[1])
					dimensionwiseorderedpolytopes.add(node)
			else:
				node = orderedarraylistnodeorig(polytope.polytopparts[0],polytope.polytopparts[1])
				dimensionwiseorderedpolytopes.add(node)
		sorted_list = list(dimensionwiseorderedpolytopes)
		sorted_list.sort(key = letter_cmp_key)
		orderedpolytoparraylist.append(sorted_list)
	return orderedpolytoparraylist

dim = 5
datasize = 70
inputpoints = tadasets.dsphere(n=datasize, d=dim-1, r=1, noise=0)  # generate d-sphere datatset
letter_cmp_key = cmp_to_key(letter_cmp)    #comparison function
#polytopaltree = createpolytopaltree(inputpoints)
polytopaltree = createnormalpolytopaltree(inputpoints)
polytopaltree = assignweights(polytopaltree)
#polytopaltree = disintegrate(polytopaltree)
for x in polytopaltree:
	print(len(x))
input()
bettieTable = set()
fastpersistance(polytopaltree)
dim0=0
dim1=0
dim2=0
dim3=0
dim4=0
for x in bettieTable:
	#print(x)
	if(x[0]==0):
		dim0= dim0+1
	if(x[0]==1):
		dim1= dim1+1
	if(x[0]==2):
		dim2= dim2+1
	if(x[0]==3):
		dim3= dim3+1
	if(x[0]==3):
		dim4= dim4+1
print("Dimemnsion 0 ::",dim0)
print("Dimemnsion 1 ::",dim1)
print("Dimemnsion 2 ::",dim2)
print("Dimemnsion 3 ::",dim3)
print("Dimemnsion 4 ::",dim4)
	
