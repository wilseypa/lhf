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


maxepsilon = 2
distancematrix = []	  #initializae distance matrix

	
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

	

def computedistancematrix(points):    #compute Distance Matrix
	global distancematrix    #initializae distance matrix
	distancematrix = []  
	for x in range(0,len(points)):
		distance = []
		for y in range(x,len(points)):
			di = math.dist(points[x],points[y])
			distance.append(di)
		distancematrix.append(distance)

def distanceindex(x,y): # compute Distance by index b/w two points
	if(x>y):
		return distancematrix[y][x-y]
	else:
		return distancematrix[x][y-x]

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
				
def minimumspanningtree(edges1):
	# using Kruskal's algorithm to find the cost of Minimum Spanning Tree
	res = 0
	pivots = []
	mstSize = 0
	ds = DisjointSet(datasize)
	for x in reversed(edges1):
		if ds.find(x.simplex[0]) != ds.find(x.simplex[1]):
			ds.union(x.simplex[0], x.simplex[1])
			res += x.weight
			mstSize +=1
			pivots.append(x)
			bettitableentry = [0,0,float(f'{x.weight:.6f}')]
			bettieTable.append(tuple(bettitableentry))
			if(mstSize >= len(datapoints)-1):
				for i in range(0,len(datapoints)):
					if(ds.find(i) == i):
						bettitableentry = [0,0,maxepsilon]
						bettieTable.append(tuple(bettitableentry))
				return pivots
	return pivots

def getDimEdges(dimension):
	edg =  Complex[dimension]
	return edg
	
def getAllCofacets(e,dimension):    #will update this function a little
	returnlist = []
	for x in Complex[dimension+1]:
		if(len(intersection(x.simplex,e.simplex)) == len(e.simplex)):
			returnlist.append(x)
	return returnlist

def intersection(lst1, lst2):  #return intersection of two lists
    return list(set(lst1) & set(lst2))
 
def union(lst1,lst2):        #return union of two lists
	return list(set().union(lst1 , lst2))


def persistenceByDimension( edges, pivots, dimension):
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
					if(e.weight < pivot.weight):
						bettitableentry = [dimension,float(f'{min(pivot.weight, e.weight):.6f}'),float(f'{max(pivot.weight, e.weight):.6f}')]
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

def fastpersistance(simplicialcomplex):
	vertices = getDimEdges(0)
	edges = getDimEdges(1)
	pivots  = minimumspanningtree(edges)
	for d  in range(1,dim):
		if(d != 1):
			 edges = getDimEdges(d)
		pivots = persistenceByDimension(edges, pivots, d)
	return

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
	
def simplexweight(simplex):
	if(len(simplex)==1):
		return 0
	facets = generatesimplexfacets(simplex,2)
	dist = 0
	for edge in facets:
		if(dist < distanceindex(edge[0],edge[1])):
			dist = distanceindex(edge[0],edge[1])
	return dist
	
def generatesimplexfacets(poly,dimen):  #generate simplices of dimension dimen
	tp = list(combinations(poly, dimen))
	listreturn = []
	for x in tp:
		listreturn.append(list(x))
	return listreturn
	
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
	
def createSimplexTree(inputpoints):
	computedistancematrix(inputpoints)
	DelaunayComplex = [set() for x in range(len(inputpoints[0])+1)]
	triangulation =  Delaunay(inputpoints).simplices
	for x in triangulation:
		for y in range(len(triangulation[0])+1):
			simplexes = itertools.combinations(x, y)
			for simplex in simplexes:
				DelaunayComplex[len(simplex)-1].add(simplex)
	DelaunayComplex = assignweights(DelaunayComplex)
	return DelaunayComplex
	

datapoints = []

datapoints = tadasets.dsphere(n=40, d = 2, r =1 , noise=0.001)	

letter_cmp_key = cmp_to_key(letter_cmp)    #comparison function

Complex = createSimplexTree(datapoints)

'''
for x in Complex:
	for y in x:
		print(y.simplex,y.weight)
'''
datasize = len(datapoints)
dim=len(datapoints[0])
	
bettieTable = []
fastpersistance(Complex)




######################################### Outline for PwPH Computaion ################################

# 1. Identify Initial epsilon Value (e0) to identify conquarbale pseudovertices.
# 2. Compute PH on each PV separtely and store the PIs.
# 3. Identify the alpha_value to identify shape boundary of each PV. 
# 4. Keep the boundary points and remove the interior points from the Alpha shape to obtain reduced point cloud P'.
# 5. Find the subsequent epsilon value so that the next set of PV's are conquarble.
# 6. Compute PH on next PV's and store the valid PIs. Need to remove the invalid PIs (Will be done by identifying the death simplex centroid inclusion in the alpha shape).
# 7. Repeate from step 3 until the Point cloud is small enough to be conquered as whole.


######################################################################################################






dimcount=[0 for i in range(0,dim)]
table = [[] for i in range(0,dim)]
for x in bettieTable:
	table[x[0]].append(x)
	dimcount[x[0]]= dimcount[x[0]]+1
for i in range(0,dim):
	print("Dimemnsion ",i," Betti Count ::",dimcount[i])
	
writebetties(table,"output.csv")
	

dimcount=[0 for i in range(0,dim)]

table = [[] for i in range(0,dim)]

for x in bettieTable:
	table[x[0]].append(x)
	dimcount[x[0]]= dimcount[x[0]]+1
for i in range(0,dim):
	print("Dimemnsion ",i," Betti Count ::",dimcount[i])

colors = ["Orange","yellow","Green",'Blue','Red','Black',"pink"]
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
print(df)

i = 0
for d in range(0,dim):
	counter = []
	[counter.append(j+i) for j in range(0,len(table[d]))]
	i=i+len(table[d])
	df1 =  df[df["Dimension"]==d]
	plt.plot([df1["Birth"], df1["Death"]], [counter, counter],color=colors[d],linestyle='solid',linewidth=1)
plt.savefig("outputPIpolytopal.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()


for d in range(0,dim):
	df1 =  df[df["Dimension"]==d]
	plt.scatter(df1["Birth"], df1["Death"],color = colors[d])
plt.axline([0, 0], [2, 2],linewidth=1,color="black")
plt.savefig("outputBCPolytopal.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()

