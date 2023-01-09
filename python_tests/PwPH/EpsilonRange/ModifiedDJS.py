import numpy as np
import csv
import math 
import tadasets
from scipy.spatial import Delaunay
import itertools
from functools import cmp_to_key
from itertools import combinations
import heapq as hq
import pandas as pd
import persim
import matplotlib.pyplot as plt
import copy
import alphashape
from descartes import PolygonPatch
from shapely.geometry import Point
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from matplotlib.collections import LineCollection

#***********************************************************

print("Enter Maximum Conquerable Size ::")
maximumconquarablesize = int(input())
maxepsilon = 999999
distancematrix = []	  #initializae distance matrix

#***********************************************************


def generatesedges(dataindexes):  #generate simplices of dimension dimen
	tp = list(combinations(dataindexes, 2))
	listreturn = []
	for x in tp:
		weight = math.dist(datapoints[x[0]],datapoints[x[1]])
		x = x + tuple([weight])
		listreturn.append(list(x))
	return sorted(listreturn, key = lambda x: x[2])

#***********************************************************

	
class orderedarraylistnode(object):
	def __init__(self,simplex,weight):
		self.simplexparts = simplex
		simplex1 = []
		if(type(simplex[0][0])!=type(np.int32(1)) and type(simplex[0][0])!=type(1)):
			[simplex1.append(x) for y in simplex for x in y[0] if x not in simplex1]
		else:
			[simplex1.append(x) for x in simplex[0] if x not in simplex1]
		simplex1 = sorted(simplex1)
		self.simplex = sorted(simplex1,reverse=False)
		self.weight = weight
	
	def __hash__(self):
		return hash((tuple(self.simplex), self.weight))
	
	def __getitem__(self, item):
		return [self.simplex,self.weight]
	
	def __lt__(self, nxt):
		if self.weight < nxt.weight:
			return True
		elif self.weight == nxt.weight:
			for (x,y) in zip(reversed(self.simplex),reversed(nxt.simplex)):
				if(x < y):
					return True
				elif(x==y):
					continue
				else:
					return False
		return False
	def __eq__(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		return self.simplex == other.simplex and self.weight == other.weight

#***********************************************************

def computedistancematrix(points):    #compute Distance Matrix
	global distancematrix    #initializae distance matrix
	distancematrix = []  
	for x in range(0,len(points)):
		distance = []
		for y in range(x,len(points)):
			di = math.dist(points[x],points[y])
			distance.append(di)
		distancematrix.append(distance)

#***********************************************************

def distanceindex(x,y): # compute Distance by index b/w two points
	if(x>y):
		return distancematrix[y][x-y]
	else:
		return distancematrix[x][y-x]

#***********************************************************

def writebetties(table,filename):
	birth = []
	death = []
	dimension = []

	for d in range(0,dim):
		for x in table[d]:
			dimension.append(x[0])
			birth.append(x[1])
			death.append(x[2])
	df = pd.DataFrame(list(zip(dimension,birth,death)),columns =['dimension','birth','death'])
	df.to_csv(filename,index=False,header=False)

#***********************************************************

class DisjointSet:
    def __init__(self, indexes,p,r):
        self.original = copy.deepcopy(indexes)
        self.parent = p #indexes
        self.rank = r# [1 for _ in indexes]
    def getPVs(self,indexes):
        djs = [[] for i in range(len(datapoints))]
        for i in indexes:
            djs[self.find(i)].append(i)
        djs = [ele for ele in djs if ele != []]
        return djs
    # make a and b part of the same component
    # union by rank optimization
    def union(self, a, b):
        pa = self.find(a)
        pb = self.find(b)
        if pa == pb: return
        if self.rank[self.original.index(pa)] > self.rank[self.original.index(pb)]:
            self.parent[self.original.index(pb)] = pa
            self.rank[self.original.index(pa)] += self.rank[self.original.index(pb)]
        else:
            self.parent[self.original.index(pa)] = pb
            self.rank[self.original.index(pb)] += self.rank[self.original.index(pa)]
    # find the representative of the 
    # path compression optimization
    def find(self, a):
        if self.parent[self.original.index(a)] == a:
            return a
        
        self.parent[self.original.index(a)] = self.find(self.parent[self.original.index(a)])
        return self.parent[self.original.index(a)]
        
    def findrank(self, a):
        return self.rank[a]
    
    def SameSet(n1, n2):
        return self.find(n1) == self.find(n2)

#***********************************************************

def merge(ds1,ds2):
	indexes  = ds1.original + ds2.original
	parents = ds1.parent + ds2.parent
	rank = ds1.rank + ds2.rank
	return DisjointSet(indexes,parents,rank)

#***********************************************************
	
def growFullMST(points):
	ds = DisjointSet(points,points,[1 for i in points])
	edges = generatesedges(points)
	mstSize = 0
	mstedges = []
	for x in edges:
		if ds.find(x[0]) != ds.find(x[1]):
			ds.union(x[0], x[1])
			mstedges.append([x[0],x[1]])
			mstSize +=1
		if(mstSize >= len(points)-1):
			return ds,mstedges
	return ds,mstedges

#***********************************************************

def growPVs(updatedVertices,ds):
	# using Kruskal's algorithm to find the cost of Minimum Spanning Tree
	res = 0
	mstSize = 0
	previousweight = 0
	mstedges = []
	edges = generatesedges(updatedVertices)
	for x in edges:
		if ds.find(x[0]) != ds.find(x[1]):
			ds.union(x[0], x[1])
			mstedges.append([x[0],x[1]])
			mstSize +=1
			maxmstsize = max([ds.findrank(i) for i in range(0,len(updatedVertices))])
			bettitableentry = [0,0,float(f'{x[2]:.6f}'),tuple([x[0],x[1]])]
			bettieTable.append(tuple(bettitableentry))
			if(mstSize >= len(updatedVertices)-1 or maxmstsize > maximumconquarablesize):
				return ds,x[2],mstedges
	return ds, 9999,mstedges
\
#***********************************************************
 
def circumRadius(simplex,points):
	matA = [[] for i in range(len(simplex))]
	matACap = [[] for i in range(len(simplex)+1)]
	ii=0
	for i in simplex:
		matACap[ii+1].append(1)
		for j in simplex:
			if(math.dist(points[i],points[j])!=0):
				matA[ii].append(pow(math.dist(points[i],points[j]),2))
				matACap[ii+1].append(pow(math.dist(points[j],points[i]),2))
			else:
				matA[ii].append(pow(math.dist(points[j],points[i]),2))
				matACap[ii+1].append(pow(math.dist(points[i],points[j]),2))
		ii= ii+1
	matACap[0].append(0)
	for i in simplex:
		matACap[0].append(1)
	return math.sqrt(-(np.linalg.det(np.array(matA))/(2*np.linalg.det(np.array(matACap)))))

#***********************************************************

def intersection(lst1, lst2):  #return intersection of two lists
    return list(set(lst1) & set(lst2))

#***********************************************************    
    
def generatesimplexfacets(simplex,dimen):  #generate simplices of dimension dimen
	tp = list(combinations(simplex, dimen))
	listreturn = []
	for x in tp:
		listreturn.append(list(x))
	return listreturn
	
#***********************************************************	

def getDimEdges(dimension):
	edg =  Complex[dimension]
	return edg
	
#***********************************************************	

def getAllCofacets(e,dimension):    #will update this function a little
	returnlist = []
	for x in Complex[dimension+1]:
		if(len(intersection(x.simplex,e.simplex)) == len(e.simplex)):
			returnlist.append(x)
	return returnlist

#***********************************************************	

def alphaBoundary(alphaShape,isolatedVertices):
	boundaryVertices = isolatedVertices
	alphahull = []
	for x in alphaShape:
		facets = generatesimplexfacets(x,len(x)-1)
		for f in facets:
			valid = True
			for y in alphaShape:
				if(valid==False):
					break
				if(x !=y):
					if(len(intersection(f,y))==len(f)):
						valid=False
						break
			if(valid==True):
				if(sorted(f) not in alphahull):
					alphahull.append(sorted(f))
	for x in alphahull:
		for y in x:
			if y not in boundaryVertices:
				boundaryVertices.append(y)
	return boundaryVertices

#***********************************************************

def centroid(simplex,newpoints):
	cntr = []
	for i in range(len(newpoints[0])):
		coord = 0
		for x in simplex:
			coord = coord + newpoints[x][i]
		coord = coord/len(simplex)
		cntr.append(coord)
	return cntr

#***********************************************************

def persistenceByDimension( edges, pivots, dimension,newpoints):
	pivotcounter = len(pivots)-1
	edges.sort(key = letter_cmp_key)
	nextPivots = []	
	v = {} 
	pivotPairs = {}
	for e in edges:
		if(pivotcounter < 0 or pivots[pivotcounter].simplex != e.simplex):
			faceList = getAllCofacets(e,dimension)
			columnV = []
			columnV.append(e)
			hq.heapify(faceList)
			i = 1
			while(True):
				hq.heapify(faceList)
				if(i%10000==0):
					print(i)
				i=i+1
				while(faceList!=[]):
					hq.heapify(faceList)
					pivot =	hq.heappop(faceList)
					if(faceList!=[] and pivot.simplex==faceList[0].simplex):
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
						centroidsimplex = centroid(pivot.simplex,newpoints)
						bettitableentry = [dimension,float(f'{min(pivot.weight, e.weight):.6f}'),float(f'{max(pivot.weight, e.weight):.6f}'),tuple(centroidsimplex)]
						bettieTable.append(tuple(bettitableentry))
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

#***********************************************************
countPHPVS = 0
def fastpersistance(simplicialcomplex,newpoints):
	global countPHPVS
	countPHPVS = countPHPVS + 1
	print("Computing PH of PV", countPHPVS, len(newpoints))

	vertices = getDimEdges(0)
	edges = getDimEdges(1)
	pivots  = minimumspanningtree(edges,newpoints)
	for d  in range(1,dim):
		if(d != 1):
			 edges = getDimEdges(d)
		pivots = persistenceByDimension(edges, pivots, d,newpoints)
	return

#***********************************************************

def minimumspanningtree(edges1,datapoints):
	# using Kruskal's algorithm to find the cost of Minimum Spanning Tree
	res = 0
	pivots = []
	mstSize = 0
	ds = DisjointSet([i for i in range(datasize)],[i for i in range(datasize)],[1 for i in range(datasize)])
	for x in reversed(edges1):
		if ds.find(x.simplex[0]) != ds.find(x.simplex[1]):
			ds.union(x.simplex[0], x.simplex[1])
			res += x.weight
			mstSize +=1
			pivots.append(x)
			#bettitableentry = [0,0,float(f'{x.weight:.6f}'),tuple([1000,1000])]
			#bettieTable.append(tuple(bettitableentry))
			if(mstSize >= len(datapoints)-1):
				#for i in range(0,len(datapoints)):
					#if(ds.find(i) == i):
						#bettitableentry = [0,0,maxepsilon,tuple([1000,1000])]
						#bettieTable.append(tuple(bettitableentry))
				return pivots
	return pivots

#***********************************************************

def letter_cmp(a, b):
    if a.weight > b.weight:
        return -1
    elif a.weight == b.weight:
        if a.simplex > b.simplex:
            return 1
        else:
            return -1
    else:
        return 1 
        
#***********************************************************
	
def pointInsideSimplex(simplex,point,points):
	matT = []
	matPv = []
	for x  in range(len(simplex)-1):
		tempmat = []
		for j in range(len(points[0])):
			tempmat.append(points[x][j]-points[len(simplex)-1][j])
		matT.append(tempmat)
	for j in range(len(points[0])):
		tempmat = []
		tempmat.append(point[j]-points[len(simplex)-1][j])
		matPv.append(tempmat)
	transposematT = [[0 for j in range(len(matT[0]))] for i in range(len(matT))]
	for i in range(len(matT)):
		for j in range(len(matT[0])):
			transposematT[j][i]= matT[i][j]
	lam = np.matmul(np.linalg.inv(transposematT),matPv)
	outside = False;
	sum1 = 0;
	for x in lam:
		sum1 = sum1 + x[0]
		if(x[0]<0):
			outside = True
			break
	if(sum1 > 1):
		outside = True
	
	return not outside
 
#***********************************************************

def simplexweight(simplex):
	if(len(simplex)==1):
		return 0
	facets = generatesimplexfacets(simplex,2)
	dist = 0
	for edge in facets:
		if(dist < distanceindex(edge[0],edge[1])):
			dist = distanceindex(edge[0],edge[1])
	return dist

#***********************************************************

def assignweights(simplextree):
	orderedsimplexarraylist = []
	for dimensionlist in simplextree:
		dimensionwiseorderedpolytopes = set()
		for simplex in dimensionlist:	
			maxedge = 0
			smplex = []
			smplex.append([simplex,simplexweight(simplex)])
			node = orderedarraylistnode(smplex,simplexweight(simplex))
			dimensionwiseorderedpolytopes.add(node)
		sorted_list = list(dimensionwiseorderedpolytopes)
		sorted_list.sort(key = letter_cmp_key)
		orderedsimplexarraylist.append(sorted_list)
	return orderedsimplexarraylist

#***********************************************************

def createSimplexTree(inputpoints,epsilon):
	DelaunayComplex = [set() for x in range(len(inputpoints[0])+1)]
	triangulation =  Delaunay(inputpoints).simplices
	for x in triangulation:
		for y in range(len(inputpoints[0])+1):
			y= y+1
			simplexes = itertools.combinations(x, y)
			for simplex in simplexes:
				if(simplexweight(simplex)<=(epsilon)):
					DelaunayComplex[len(simplex)-1].add(simplex)
	DelaunayComplex = assignweights(DelaunayComplex)
	return DelaunayComplex

#***********************************************************

def createVRComplexTree(inputpoints,epsilon):
	VRComplex = [set() for x in range(len(inputpoints[0])+1)]
	triangulation = range(len(inputpoints))
	for y in range(len(inputpoints[0])+1):
		y = y+1
		simplexes = itertools.combinations(triangulation, y)
		for simplex in simplexes:
			if(simplexweight(simplex)<=(epsilon+.0000001)):
				VRComplex[len(simplex)-1].add(simplex)
	VRComplex = assignweights(VRComplex)
	return VRComplex

#***********************************************************

def computeAlphaShape(vertices,afv):
	isolatedVertices = []
	uniqueIndices = set()
	maxedge = 0
	points = [datapoints[i] for i in vertices]
	triangulation = Delaunay(points).simplices
	alphaShape = []
	for simp in triangulation:
		if(circumRadius(simp,points) < afv):
			simplex = []
			for y in simp:
				simplex.append(vertices[y])
				uniqueIndices.add(vertices[y])
			alphaShape.append(simplex)
			weight = simplexweight(simplex)
			if(maxedge<weight):
				maxedge = weight
	for x in vertices:
		if x not in uniqueIndices:
			isolatedVertices.append(x)
	return alphaShape,isolatedVertices,maxedge

#***********************************************************

def filterPIs(bettieTable,alphashapesimplices,epsilon):
	prunedPis = []
	for x in bettieTable:
		inside = False
		for y in alphashapesimplices:
			if pointInsideSimplex(y,x[3],datapoints):
				inside = True
				break
		if not inside and x[2]>epsilon:
			prunedPis.append(x)
	return prunedPis

#***********************************************************
letter_cmp_key = cmp_to_key(letter_cmp)    #comparison function

#datapoints = tadasets.dsphere(n=200, d = 1, r =4 , noise=0.1)	
#datapoints2 = tadasets.dsphere(n=200, d = 1, r =2 , noise=0.1)

#datapoints = np.vstack((datapoints,datapoints2))
#filename = "PathBased"
#filename = "ZahnsCompund"
#filename = "D31"
print("Enter File Name")
filename = input()
datapoints = np.loadtxt(filename+".csv",delimiter=",", dtype=float)

print("Enter Complex to Use:: (Del/VR)")
complextype = input()
#datapoints = np.loadtxt("ZahnsCompund.csv",delimiter=",", dtype=float)
#datapoints = np.loadtxt("D31.csv",delimiter=",", dtype=float)

dim = len(datapoints[0])
datasize = len(datapoints)
#******************************************
updatedPVs = [[i] for i in range(len(datapoints))]
finalfigure = []
globalseparatebettieTable = []
globalbettieTable = set()
alphashapesimplices = set()
alpha_shape = []
computedistancematrix(datapoints)
while True:
	print("********")
	reduceddatasize = 0
	updatedVertices = []
	updatedVertices += updatedPVs[0]
	globalmstedges = []
	reduceddatasize += len(updatedPVs[0])
	globalDS,mstedge = growFullMST(updatedPVs[0])
	if(mstedge):
		globalmstedges.append(mstedge)
	for i in range(1,len(updatedPVs)):
		updatedVertices += updatedPVs[i]
		ds,mstedge = growFullMST(updatedPVs[i])
		globalDS = merge(globalDS,ds)
		reduceddatasize += len(updatedPVs[i])
		if(mstedge):
			globalmstedges.append(mstedge)
		'''j = 0
		for x,y,z in zip(ds.parent,ds.rank,ds.original):
			print(j," Parent ::", x," Rank ::", y, " Indices :: ", z)
			j += 1
		print("********")
		'''
	getbetties = set()
	print("Size::",reduceddatasize)
	bettieTable = []
	pseudoVertices = growPVs(updatedVertices,globalDS)
	for x in bettieTable:
		globalbettieTable.add(x)
		getbetties.add(x)

	'''
	i = 0
	for x,y,z in zip(pseudoVertices[0].parent,pseudoVertices[0].rank,pseudoVertices[0].original):
		print(i," Parent ::", x," Rank ::", y, " Indices :: ", z)
		i += 1
	print("********")
	print("Vertices ", updatedVertices)
	print("PseudoVertices  ")
	'''
	list2= []
	newmstedges = pseudoVertices[2]
	PVs = pseudoVertices[0].getPVs(updatedVertices)
	print("Total PVs ", len(PVs))
	alpha_value = pseudoVertices[1]
	print("Alpha Value  ",alpha_value)
	if len(PVs)==1:
		pv = PVs[0]
		for x in alpha_shape:
			alphashapesimplices.add(tuple(x))
		#****************************************************
		alpha_shape,isolatedVertices,PHepsilon = computeAlphaShape(PVs[0],maxepsilon)
		list2.append([alpha_shape,isolatedVertices])
		finalfigure.append([list2,globalmstedges,newmstedges])
		#****************************************************
		
		bettieTable = []
		PVpoints = [datapoints[i] for i in pv]
		if complextype == "Del" or complextype == "del" or complextype == "Delaunay" or complextype == "delaunay":
			Complex = createSimplexTree(PVpoints,PHepsilon)
		else:
			Complex = createVRComplexTree(PVpoints,maxepsilon)
		fastpersistance(Complex,PVpoints)
		
		prunedPIs = filterPIs(bettieTable,alphashapesimplices,alpha_value)
		for x in prunedPIs:
			globalbettieTable.add(x)
			getbetties.add(x)
		globalseparatebettieTable.append(getbetties)
		break
	updatedPVs = []
	for pv in PVs:
		updatedPV = []
		if (len(pv)>(dim+1)):
			for x in alpha_shape:
				alphashapesimplices.add(tuple(x))
			#****************************************************
			alpha_shape,isolatedVertices,PHepsilon = computeAlphaShape(pv,alpha_value)
			boundaryVertices = alphaBoundary(alpha_shape,isolatedVertices)
			list2.append([alpha_shape,isolatedVertices])
		    #****************************************************
			for pt in boundaryVertices:
				updatedPV.append(pt)
			
			bettieTable = []
			PVpoints = [datapoints[i] for i in pv]
			if complextype == "Del" or complextype == "del" or complextype == "Delaunay" or complextype == "delaunay":
				Complex = createSimplexTree(PVpoints,PHepsilon)
			else:
				Complex = createVRComplexTree(PVpoints,PHepsilon)
			fastpersistance(Complex,PVpoints)
			prunedPIs = filterPIs(bettieTable,alphashapesimplices,alpha_value)
			for x in prunedPIs:
				globalbettieTable.add(x)
				getbetties.add(x)
		else:
			#5 Collect Isolated point for reduced dataset
			for pt in pv:
				updatedPV.append(pt)
			list2.append([[],updatedPV])
		updatedPVs.append(updatedPV)
	globalseparatebettieTable.append(getbetties)
	finalfigure.append([list2,globalmstedges,newmstedges])

#**************************************
print("Done Post Processing")

'''
number1 = int(math.sqrt(len(finalfigure)))
number2 = int(((len(finalfigure))/number1)+1)
rows, cols = number1,number2
fig, ax = plt.subplots(rows, cols,sharex='col', sharey='row')

i = -1
for xx in range(rows):
	for yy in range(cols):
		if(i==-1):
			x,y = zip(*datapoints)
			ax[xx,yy].scatter(x,y,s=0.5,color = "blue")
			i = i+1
		else:
			plotdata = finalfigure[i]
			plot1 = plotdata[0]
			mst2 = np.array(plotdata[2])
			mst1 = plotdata[1]
			if(mst1):
				for y in mst1:
					edg = np.array(y)
					lc1 = LineCollection(datapoints[edg],linewidths=(0.5))
					ax[xx,yy].add_collection(lc1)
			if(mst2.size>0):
				lc2 = LineCollection(datapoints[mst2],linewidths=(0.5))
				ax[xx,yy].add_collection(lc2)
			for x in plot1:
				alphashape = x[0]
				isolatedVertices = x[1]
				for j in alphashape:
					X = np.array([list(datapoints[i]) for i in j])
					ax[xx,yy].scatter(X[:, 0], X[:, 1],s=0.5, color = "blue")
					t1 = plt.Polygon(X[:3,:], color="blue",alpha = 0.2,lw = 0.2)
					ax[xx,yy].add_patch(t1)
				IV = np.array([list(datapoints[i]) for i in isolatedVertices])
				if(IV.size > 0):
					ax[xx,yy].scatter(IV[:, 0], IV[:, 1],s=0.5, color = "red")
			i = i + 1
		if(i >= len(finalfigure)):
			break
	if(i >= len(finalfigure)):
		break
plt.savefig("Subsequent"+filename+".pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()


datapoints1 = tadasets.dsphere(n=10, d = 1, r =4 , noise=0.3)	
datapoints2 = tadasets.dsphere(n=10, d = 1, r =2 , noise=0.3)

datapoints = np.vstack((datapoints1,datapoints2))
datasize = len(datapoints)
dataindexes = [i for i in range(datasize)]
even = [i for i in range(datasize) if i<5]
odd = [i for i in range(datasize) if i<10 and i>4]
odd1 = [i for i in range(datasize) if i<20 and i>9]

datapointso = [datapoints[i] for i in odd]
datapointse = [datapoints[i] for i in even]
datapointso1 = [datapoints[i] for i in odd1]

x,y = zip(*datapoints1)
plt.scatter(x,y,s=10,color = "blue")
x,y = zip(*datapoints2)
plt.scatter(x,y,s=10,color = "green")
plt.show()	

ds1 = DisjointSet(odd,odd,[1 for i in odd])

edges1 = generatesedges(odd)

for x in edges1:
	if ds1.find(x[0]) != ds1.find(x[1]):
		ds1.union(x[0], x[1])
for x,y,z in zip(ds1.parent,ds1.rank,ds1.original):
	print(" Parent ::", x," Rank ::", y, " Indices :: ", z)
print("********")

ds2 = DisjointSet(even,even,[1 for i in even])

edges1 = generatesedges(even)

for x in edges1:
	if ds2.find(x[0]) != ds2.find(x[1]):
		ds2.union(x[0], x[1])
for x,y,z in zip(ds2.parent,ds2.rank,ds2.original):
	print(" Parent ::", x," Rank ::", y, " Indices :: ", z)
print("********")

ds3 = DisjointSet(odd1,odd1,[1 for i in odd1])

edges1 = generatesedges(odd1)

for x in edges1:
	if ds3.find(x[0]) != ds3.find(x[1]):
		ds3.union(x[0], x[1])
for x,y,z in zip(ds3.parent,ds3.rank,ds3.original):
	print(" Parent ::", x," Rank ::", y, " Indices :: ", z)
print("********")


mergeddiskointset = merge(ds1,ds2)
mergeddiskointset  = merge(mergeddiskointset ,ds3)

edges1 = generatesedges(dataindexes)

for x in edges1:
	print(x)
	if mergeddiskointset.find(x[0]) != mergeddiskointset.find(x[1]):
		mergeddiskointset.union(x[0], x[1])
	
for x,y,z in zip(mergeddiskointset.parent,mergeddiskointset.rank,mergeddiskointset.original):
	print(" Parent ::", x," Rank ::", y, " Indices :: ", z)
print("********")
'''
#globalseparatebettieTable.append(globalbettieTable)
#print("Length of Iteration",len(globalseparatebettieTable))
number1 = int(math.sqrt(len(globalseparatebettieTable)))
number2 = int(((len(globalseparatebettieTable))/number1)+1)
rows, cols = number1,number2
fig, ax = plt.subplots(rows,cols,sharex='col', sharey='row')
fig1, ax1 = plt.subplots(rows,cols,sharex='col', sharey='row')
iteration=0
maxvalue = 0
for x in distancematrix:
	if maxvalue < max(x):
		maxvalue = max(x)
bettieTable = []
for globalbettieTable in globalseparatebettieTable:
	xx = int(iteration/cols)
	yy = int(iteration%cols)
	for it in globalbettieTable:
		bettieTable.append(it)
	dimcount=[0 for i in range(0,dim)]
	table = [[] for i in range(0,dim)]
	for x in bettieTable:
		if(x[2] != maxepsilon):
			table[x[0]].append(x)
			dimcount[x[0]]= dimcount[x[0]]+1
	
	writebetties(table,"output"+str(i)+filename+".csv")
	

	dimcount=[0 for i in range(0,dim)]

	table = [[] for i in range(0,dim)]

	for x in bettieTable:
		if(x[2] != maxepsilon):
			table[x[0]].append(x)
			dimcount[x[0]]= dimcount[x[0]]+1
	#for i in range(0,dim):	
	#	print("Dimemnsion ",i," Betti Count ::",dimcount[i])

	colors = ["Green",'Blue','Red','Black',"Orange","yellow","pink"]
	Dimension = [x[0] for x in table[0]]
	Birth = [x[1] for x in table[0]]
	Death = [x[2] for x in table[0]]
	df = pd.DataFrame(list(zip(Dimension,Birth,Death)),columns =['Dimension','Birth','Death'])

	for d in range(1,dim):
		Dim = [x[0] for x in table[d]]
		B = [x[1] for x in table[d]]
		D = [x[2] for x in table[d]]
		dftemp = pd.DataFrame(list(zip(Dim,B,D)),columns =['Dimension','Birth','Death'])
		df = pd.concat([df, dftemp], axis=0)

	df = df.sort_values(by=['Death'])

	i = 0
	for d in range(0,dim):
		counter = []
		[counter.append(j+i) for j in range(0,len(table[d]))]
		i=i+len(table[d])
		df1 =  df[df["Dimension"]==d]
		ax[xx][yy].plot([df1["Birth"], df1["Death"]], [counter, counter],color=colors[d],linestyle='solid',linewidth=1)
		ax[xx][yy].set_xlim([0,maxvalue])
	for d in range(0,dim):
		df1 =  df[df["Dimension"]==d]
		ax1[xx][yy].scatter(df1["Birth"], df1["Death"],color = colors[d],s=0.5)	
		ax1[xx][yy].axline([0, 0], [maxvalue, maxvalue],linewidth=1,color="black")
		ax1[xx][yy].set_xlim([0,maxvalue])

	iteration += 1
fig.savefig("outputBarcode"+filename+".pdf", bbox_inches = 'tight',pad_inches = 0)
fig.show()
	
fig1.savefig("outputPI"+filename+".pdf", bbox_inches = 'tight',pad_inches = 0)
fig1.show()
