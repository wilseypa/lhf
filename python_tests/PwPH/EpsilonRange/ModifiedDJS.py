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

def generatesedges(dataindexes):  #generate simplices of dimension dimen
	tp = list(combinations(dataindexes, 2))
	listreturn = []
	for x in tp:
		weight = math.dist(datapoints[x[0]],datapoints[x[1]])
		x = x + tuple([weight])
		listreturn.append(list(x))
	return sorted(listreturn, key = lambda x: x[2])
	
class DisjointSet:
    def __init__(self, indexes,p,r):
        self.original = copy.deepcopy(indexes)
        self.parent = p #indexes
        self.rank = r# [1 for _ in indexes]
    def getPVs(self,indexes):
        djs = [[] for i in indexes]
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
	
datapoints1 = tadasets.dsphere(n=10, d = 1, r =4 , noise=0.3)	
datapoints2 = tadasets.dsphere(n=10, d = 1, r =2 , noise=0.3)

datapoints = np.vstack((datapoints1,datapoints2))
datasize = len(datapoints)
dataindexes = [i for i in range(datasize)]
even = [i for i in range(datasize) if i%2==0]
odd = [i for i in range(datasize) if i%2==1]

datapointso = [datapoints[i] for i in odd]
datapointse = [datapoints[i] for i in even]

x,y = zip(*datapoints1)
plt.scatter(x,y,s=10,color = "blue")
x,y = zip(*datapoints2)
plt.scatter(x,y,s=10,color = "green")
plt.show()	

ds1 = DisjointSet(odd,odd,[1 for i in odd])

edges1 = generatesedges(odd)

for x in edges1:
	print(x)
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

mergeddiskointset = merge(ds1,ds2)
for x,y,z in zip(mergeddiskointset.parent,mergeddiskointset.rank,mergeddiskointset.original):
	print(" Parent ::", x," Rank ::", y, " Indices :: ", z)
print("********")

