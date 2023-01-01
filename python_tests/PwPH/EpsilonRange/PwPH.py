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



maxepsilon = 9999
distancematrix = []	  #initializae distance matrix
maximumconquarablesize = 100
	
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
    def getPVs(self,n):
        djs = [[] for i in range(n)]
        for i in range(n):
            djs[self.find(i)].append(i)
        djs = [ele for ele in djs if ele != []]
        return djs
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
        
    def findrank(self, a):
        return self.rank[a]
    
    def SameSet(n1, n2):
        return self.find(n1) == self.find(n2)


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
				
def minimumspanningtree(edges1,datapoints):
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
	
def minimumspanningtreeforinitialization(edges1,datapoints):
	# using Kruskal's algorithm to find the cost of Minimum Spanning Tree
	res = 0
	mstSize = 0
	ds = DisjointSet(len(datapoints))
	previousweight = 0
	for x in reversed(edges1):
		if ds.find(x.simplex[0]) != ds.find(x.simplex[1]):
			previousds = copy.deepcopy(ds)
			ds.union(x.simplex[0], x.simplex[1])
			res += x.weight
			mstSize +=1
			maxmstsize = max([ds.findrank(i) for i in range(0,len(datapoints))])
			if(mstSize >= len(datapoints)-1 or maxmstsize > maximumconquarablesize):
				#return previousds,previousweight,x.weight
				return ds,x.weight,previousweight
			previousweight = x.weight
	return ds
	

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

def fastpersistance(simplicialcomplex,newpoints):
	vertices = getDimEdges(0)
	edges = getDimEdges(1)
	pivots  = minimumspanningtree(edges,newpoints)
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
	
def pointInsideSimplex(simplex,point,points):
	i = 0
	matT = []
	matPv = []
	print(simplex)
	print(point)
	print(points)
	for x  in simplex:
		if i!=1:
			tempmat = []
			for j in range(len(points[0])):
				tempmat.append(points[x][1]-points[x][j])
			matT.append(tempmat)
		if(i==1):
			for j in range(len(points[0])):
				tempmat = []
				tempmat.append(point[j]-points[x][j])
				matPv.append(tempmat)
		i = i + 1
	print(matT)
	transposematT = [[0 for j in range(len(matT[0]))] for i in range(len(matT))]
	for i in range(len(matT)):
		for j in range(len(matT[0])):
			transposematT[j][i]= matT[i][j]
	print(transposematT)
	lam = np.matmul(transposematT,matPv)
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
 
def computeAlphaShape(points,afv):
	triangulation = Delaunay(points).simplices
	alphaShape = []
	for simp in triangulation:
		if(circumRadius(simp,points) < afv):
			alphaShape.append(simp)
	return alphaShape
		

def alphaBoundary(alphaShape,isolatedVertices):
	alphahull= []
	boundaryVertices = isolatedVertices
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
	return alphahull,boundaryVertices
		
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

def createVRComplexTree(inputpoints,epsilon,dimension):
	computedistancematrix(inputpoints)
	VRComplex = [set() for x in range(len(inputpoints[0])+1)]
	triangulation = range(len(inputpoints))
	for y in range(dimension+1):
		simplexes = itertools.combinations(triangulation, y)
		for simplex in simplexes:
			VRComplex[len(simplex)-1].add(simplex)
	VRComplex = assignweights(VRComplex)
	return VRComplex

	
def alphaShapeMapAppendIndices(alpha_shape,pv,hull):
	mappedindices = []
	hullindices = []
	uniqueIndices = set()
	isolatedVertices = []
	for x in hull:
		simplex = []
		for y in x:
			simplex.append(pv[y])
		hullindices.append(simplex)
	
	for x in alpha_shape:
		simplex = []
		for y in x:
			simplex.append(pv[y])
			uniqueIndices.add(pv[y])
		mappedindices.append(simplex)
	for x in pv:
		if x not in uniqueIndices:
			isolatedVertices.append(x)
	return mappedindices,isolatedVertices,hullindices

def MapPVs(Pvs,nps):
	mappedindices = []
	for x in Pvs:
		pseudov = []
		for y in x:
			pseudov.append(nps[y])
		mappedindices.append(pseudov)
	return mappedindices
	

def plotalphaShape(list1,newpoints):
	x,y = zip(*newpoints)
	plt.scatter(x,y,s=5,color = "red")
	for x in list1:
		alphashape = x[0]
		points = x[1]
		isolatedVertices = x[2]
		b = x[3]
		bv = x[4]
		hull = x[5]
		for j in alphashape:
			X = np.array([list(points[i]) for i in j])
			plt.scatter(X[:, 0], X[:, 1], s = 5, color = "blue")
			t1 = plt.Polygon(X[:3,:], color="blue",alpha = 0.3)
			plt.gca().add_patch(t1)
		for j in b:
			X = np.array([list(points[i]) for i in j])
			plt.scatter(X[:, 0], X[:, 1], s = 5, color = "green")
			t1 = plt.Polygon(X[:3,:], color="green",alpha = 1,lw = 2)
			plt.gca().add_patch(t1)
		for j in hull:
			X = np.array([list(points[i]) for i in j])
			plt.scatter(X[:, 0], X[:, 1], s = 5, color = "green")
			t1 = plt.Polygon(X[:3,:], color="pink",alpha = 1,lw = 2)
			plt.gca().add_patch(t1)
		IV = np.array([list(points[i]) for i in isolatedVertices])
		if(IV.size > 0):
			plt.scatter(IV[:, 0], IV[:, 1], s = 10, color = "red")
		bvi = np.array([list(points[i]) for i in bv])
		if(bvi.size > 0):
			plt.scatter(bvi[:, 0], bvi[:, 1], s = 20, color = "black")
	plt.show()
#print(pointInsideSimplex([0,1,2],[2,2],[[1,1],[3,1],[2,3]]))
#input()
#datapoints = []

#datapoints = tadasets.dsphere(n=50, d = 1, r =4 , noise=0.1)	
#datapoints2 = tadasets.dsphere(n=50, d = 1, r =5 , noise=0.1)

#datapoint = np.vstack((datapoints,datapoints2))

#x,y = zip(*datapoint)
#plt.scatter(x,y)
#plt.show()

#alpha_shape = computeAlphaShape(datapoint,3)
#alpha_shape = mapindices(alpha_shape,pv)
#plotalphaShape(alpha_shape,datapoint)
			
#input()

letter_cmp_key = cmp_to_key(letter_cmp)    #comparison function


'''
for x in Complex:
	for y in x:
		print(y.simplex,y.weight)
'''

	
datapoints = tadasets.dsphere(n=200, d = 1, r =4 , noise=0.3)	
datapoints2 = tadasets.dsphere(n=200, d = 1, r =2 , noise=0.3)

datapoints = np.vstack((datapoints,datapoints2))
######################################### Outline for PwPH Computaion ################################

# 1. Identify Initial epsilon Value (e0) to identify conquarbale pseudovertices.
# 2. Compute PH on each PV separtely and store the PIs.
# 3. Identify the alpha_value to identify shape boundary of each PV. 
# 4. Keep the boundary points and remove the interior points from the Alpha shape to obtain reduced point cloud P'.
# 5. Find the subsequent epsilon value so that the next set of PV's are conquarble.
# 6. Compute PH on next PV's and store the valid PIs. Need to remove the invalid PIs (Will be done by identifying the death simplex centroid inclusion in the alpha shape).
# 7. Repeate from step 3 until the Point cloud is small enough to be conquered as whole.


######################################################################################################
datasize = len(datapoints)
dim = len(datapoints[0])
newpoints = [i for i in range(datasize)]
globalbettieTable = []
finalfigure = []
while(True):
	print("Reduced Point Cloud " ,len(newpoints))
	#print(sorted(newpoints), "  ", len(newpoints))
	#for t in range(len(datapoints)):
	#	if t not in newpoints:
	#		print(t)
	nps = newpoints
	newpoints1 = [datapoints[i] for i in newpoints]
	newpoints = [datapoints[i] for i in newpoints]
	datasize = len(newpoints)
	Complex = createVRComplexTree(newpoints,maxepsilon,dim)
	edges = getDimEdges(1)
	#print(edges)
	# 1 Identifying pseduoVertices and value of cutoff Epsilon
	pseudoVertices  = minimumspanningtreeforinitialization(edges,newpoints)
	PVs = pseudoVertices[0].getPVs(datasize)
	#print(PVs)
	#print(nps)
	PVs = MapPVs(PVs,nps)
	#print(PVs)
	alpha_value = pseudoVertices[1]
	newpoints = []
	print("Subsequent Epsilons", alpha_value)
	print("Total number of Pseudo Vertices ", len(PVs), " ")
	if len(PVs)==1:
		bettieTable = []
		PVpoints = [datapoints[i] for i in PVs[0]]
		Complex = createVRComplexTree(PVpoints,maxepsilon,dim)
		fastpersistance(Complex,PVpoints)
		globalbettieTable.append(bettieTable)
		alpha_shape = computeAlphaShape(PVpoints,alpha_value)
		hull = ConvexHull(PVpoints,qhull_options="QJ").simplices
		alpha_shape,isolatedVertices,hull = alphaShapeMapAppendIndices(alpha_shape,PVs[0],hull)
		boundary, boundaryVertices = alphaBoundary(alpha_shape,isolatedVertices)
		list2.append([alpha_shape,datapoints,isolatedVertices,boundary,boundaryVertices,hull])
		finalfigure.append([list2,newpoints1])	
		break;
	list2= []
	for pv in PVs:
		if (len(pv)>dim):
			# 2 Compute PH on PV and store the results
			bettieTable = []
			PVpoints = [datapoints[i] for i in pv]
			Complex = createVRComplexTree(PVpoints,maxepsilon,dim)
			fastpersistance(Complex,PVpoints)
			# 2a Store all of the if First level PV's
			# 2b Prune Invalid PVs after first level PV's
			globalbettieTable.append(bettieTable)
		
			# 3 Compute Alpha Shape on Each PV and compute its boundary
			alpha_shape = computeAlphaShape(PVpoints,alpha_value)
			hull = ConvexHull(PVpoints,qhull_options="QJ").simplices
			alpha_shape,isolatedVertices,hull = alphaShapeMapAppendIndices(alpha_shape,pv,hull)
			boundary, boundaryVertices = alphaBoundary(alpha_shape,isolatedVertices)
			list2.append([alpha_shape,datapoints,isolatedVertices,boundary,boundaryVertices,hull])
			# 4 Consider Only boundary Vertices for next iteration
			# (Remove Interior Points)
			
			#5 Collect point for reduced dataset
			for pt in boundaryVertices:
				newpoints.append(pt)
			'''
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
			'''
		else:
			#5 Collect Isolated point for reduced dataset
			for pt in pv:
				newpoints.append(pt)
	#plotalphaShape(list2,newpoints1)
	finalfigure.append([list2,newpoints1])
	plt.show()

# Plot Final Figure

#plt.figure(figsize=(6, 4))

number = int(math.sqrt(len(finalfigure)))
rows, cols = number,number+1
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
			list1 = plot1[0]
			newpoints = plot1[1]
			x,y = zip(*newpoints)
			ax[xx,yy].scatter(x,y,s=5,color = "red")
			for x in list1:
				alphashape = x[0]
				points = x[1]
				isolatedVertices = x[2]
				b = x[3]
				bv = x[4]
				hull = x[5]
				for j in alphashape:
					X = np.array([list(points[i]) for i in j])
					ax[xx,yy].scatter(X[:, 0], X[:, 1], s = 5, color = "blue")
					t1 = plt.Polygon(X[:3,:], color="blue",alpha = 0.3)
					ax[xx,yy].add_patch(t1)
				for j in b:
					X = np.array([list(points[i]) for i in j])
					ax[xx,yy].scatter(X[:, 0], X[:, 1], s = 5, color = "green")
					t1 = plt.Polygon(X[:3,:], color="green",alpha = 1,lw = 2)
					ax[xx,yy].add_patch(t1)
				for j in hull:
					X = np.array([list(points[i]) for i in j])
					ax[xx,yy].scatter(X[:, 0], X[:, 1], s = 5, color = "green")
					t1 = plt.Polygon(X[:3,:], color="pink",alpha = 1,lw = 2)
					ax[xx,yy].add_patch(t1)
				IV = np.array([list(points[i]) for i in isolatedVertices])
				if(IV.size > 0):
					ax[xx,yy].scatter(IV[:, 0], IV[:, 1], s = 5, color = "red")
				bvi = np.array([list(points[i]) for i in bv])
				if(bvi.size > 0):
					ax[xx,yy].scatter(bvi[:, 0], bvi[:, 1], s = 5, color = "black")
			i = i + 1
			if(i >= len(finalfigure)):
				break
plt.savefig("SubsequentPVs.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()
