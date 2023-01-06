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

maximumconquarablesize = 100
epsilon = 9999999
def generatesedges(dataindexes):  #generate simplices of dimension dimen
	tp = list(combinations(dataindexes, 2))
	listreturn = []
	for x in tp:
		weight = math.dist(datapoints[x[0]],datapoints[x[1]])
		x = x + tuple([weight])
		listreturn.append(list(x))
	return sorted(listreturn, key = lambda x: x[2])


	
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

def merge(ds1,ds2):
	indexes  = ds1.original + ds2.original
	parents = ds1.parent + ds2.parent
	rank = ds1.rank + ds2.rank
	return DisjointSet(indexes,parents,rank)
	
def growFullMST(points):
	ds = DisjointSet(points,points,[1 for i in points])
	edges = generatesedges(points)
	mstSize = 0
	for x in edges:
		if ds.find(x[0]) != ds.find(x[1]):
			ds.union(x[0], x[1])
			mstSize +=1
		if(mstSize >= len(points)-1):
			return ds
	return ds

def growPVs(updatedVertices,ds):
	# using Kruskal's algorithm to find the cost of Minimum Spanning Tree
	res = 0
	mstSize = 0
	previousweight = 0
	edges = generatesedges(updatedVertices)
	for x in edges:
		if ds.find(x[0]) != ds.find(x[1]):
			ds.union(x[0], x[1])
			mstSize +=1
			maxmstsize = max([ds.findrank(i) for i in range(0,len(updatedVertices))])
			if(mstSize >= len(updatedVertices)-1 or maxmstsize > maximumconquarablesize):
				return ds,x[2]
	return ds, 9999
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

def intersection(lst1, lst2):  #return intersection of two lists
    return list(set(lst1) & set(lst2))
    
    
def generatesimplexfacets(simplex,dimen):  #generate simplices of dimension dimen
	tp = list(combinations(simplex, dimen))
	listreturn = []
	for x in tp:
		listreturn.append(list(x))
	return listreturn
	
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
def computeAlphaShape(vertices,afv):
	isolatedVertices = []
	uniqueIndices = set()

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
	for x in vertices:
		if x not in uniqueIndices:
			isolatedVertices.append(x)
	return alphaShape,isolatedVertices

#***********************************************************
	
datapoints = tadasets.dsphere(n=200, d = 1, r =4 , noise=0.1)	
datapoints2 = tadasets.dsphere(n=200, d = 1, r =2 , noise=0.1)

datapoints = np.vstack((datapoints,datapoints2))
dim = len(datapoints[0])

#******************************************
updatedPVs = [[i] for i in range(len(datapoints))]
finalfigure = []


while True:
	print()
	print("********")
	updatedVertices = []
	updatedVertices += updatedPVs[0]
	globalDS = growFullMST(updatedPVs[0])
	for i in range(1,len(updatedPVs)):
		updatedVertices += updatedPVs[i]
		ds = growFullMST(updatedPVs[i])
		globalDS = merge(globalDS,ds)
		'''j = 0
		for x,y,z in zip(ds.parent,ds.rank,ds.original):
			print(j," Parent ::", x," Rank ::", y, " Indices :: ", z)
			j += 1
		print("********")
		'''

	pseudoVertices = growPVs(updatedVertices,globalDS)
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
	PVs = pseudoVertices[0].getPVs(updatedVertices)
	print("Total PVs ", len(PVs))
	alpha_value = pseudoVertices[1]
	print("Alpha Value  ",alpha_value)
	if len(PVs)==1:		
		alpha_shape,isolatedVertices = computeAlphaShape(PVs[0],epsilon)
		list2.append([alpha_shape,isolatedVertices])
		finalfigure.append(list2)
		break
	updatedPVs = []
	for pv in PVs:
		updatedPV = []
		if (len(pv)>dim):
		    #****************************************************
			alpha_shape,isolatedVertices = computeAlphaShape(pv,alpha_value)
			boundaryVertices = alphaBoundary(alpha_shape,isolatedVertices)
			list2.append([alpha_shape,isolatedVertices])
		    #****************************************************
			for pt in boundaryVertices:
				updatedPV.append(pt)
		else:
			#5 Collect Isolated point for reduced dataset
			for pt in pv:
				updatedPV.append(pt)
		updatedPVs.append(updatedPV)
	finalfigure.append(list2)

#**************************************

number = int(math.sqrt(len(finalfigure)))
number = number+1
rows, cols = number,number
fig, ax = plt.subplots(rows, cols,sharex='col', sharey='row')

i = -1
for xx in range(rows):
	for yy in range(cols):
		if(i==-1):
			x,y = zip(*datapoints)
			ax[xx,yy].scatter(x,y,s=5,color = "blue")
			i = i+1
		else:
			plot1 = finalfigure[i]
			for x in plot1:
				alphashape = x[0]
				isolatedVertices = x[1]
				for j in alphashape:
					X = np.array([list(datapoints[i]) for i in j])
					ax[xx,yy].scatter(X[:, 0], X[:, 1], s = 5, color = "blue")
					t1 = plt.Polygon(X[:3,:], color="blue",alpha = 0.3)
					ax[xx,yy].add_patch(t1)
				IV = np.array([list(datapoints[i]) for i in isolatedVertices])
				if(IV.size > 0):
					ax[xx,yy].scatter(IV[:, 0], IV[:, 1], s = 5, color = "red")
			i = i + 1
		if(i >= len(finalfigure)):
			break
	if(i >= len(finalfigure)):
		break
plt.savefig("SubsequentPVs.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()

'''

	
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
