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
def computedistancematrix(points):    #compute Distance Matrix
	minimumdistance = 99999
	maximumdistance = 0
	global distancematrix
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

def distanceindex(x,y): # compute Distance by index b/w two points
	computedistancematrix()
	global distancematrix 
	if(x>y):
		return distancematrix[y][x-y]
	else:
		return distancematrix[x][y-x]
		
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
	return convexpartsunion


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
	def __init__(self,polytope1,polyunmapped,weight1,polytophalyspaces,polyhalfspace,chebR1,chebXc1,farthestpoint1,steriographic_projection1,faces,convexdecomposedfaces,mappedconvexdecomposedfaces,dimension,pca_x):
		self.polytopevertices = polytope1  # list of coordinates of polytope vertices
		self.polytopunmappedvertices = polyunmapped
		self.weight = weight1  #weight of polytope
		self.polytophalyspaces = polytophalyspaces  #Hyperplane System of Equations
		self.polyhalfspace = polyhalfspace  # Create polytope object
		self.chebR = chebR1  #polytope incenter
		self.chebXc = chebXc1   #polytope inradius
		self.farthestpoint = farthestpoint1 #polytope fathest point from incenter coordinates
		self.steriographic_projection = steriographic_projection1 # list of corrdinates one dimensional less (Conformal Projection)
		self.faces = faces #list of polytopal faces
		self.convexdecomposed = convexdecomposedfaces
		self.convexdecomposedfaces =  mappedconvexdecomposedfaces  #PCA Dimensionality Reduction for lower order faces
		self.dimension = dimension # vertices in highest simplex
		self.childrens = []
		self.pca_x = pca_x

class tree:
	def createheadpolytope(self,points):
		poly = points
		weight = 0
		polytophalyspaces = []
		polyhalfspace = []
		chebR = []
		chebXc = []
		farthestpoint = 9999999
		steriographic_projection = inputpoints
		faces =  Delaunay(inputpoints).simplices
		convexdecomposedfaces = iterativeconvexization(faces,len(inputpoints[0]),inputpoints)
		dimension = len(inputpoints[0])
		pca_x = inputpoints
		return polytope(poly,poly,weight,polytophalyspaces,polyhalfspace,chebR,chebXc,farthestpoint,steriographic_projection,faces,convexdecomposedfaces,convexdecomposedfaces,dimension,pca_x)

	def createpolytope(self,points,points1,pca_x1,farthest,dimension):
		poly = points
		polyunmapped = points1
		#print(points)
		polytopal_points = []

		if(farthest ==-1 or len(points1) == dimension+1):
			mappedconvexfaces = generatesimplexfacets(poly, dimension)
			convexfaces = generatesimplexfacets([i for i in range(0,len(points))], dimension)
			return polytope(poly,points1,0,[],[],[],[],-1,[],[],convexfaces,mappedconvexfaces,dimension-1,[])
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
		if(len(pca_x[0]) <= 1):
			faces=[]
			pcx = []
			for x in pca_x:
				for k in x:
					pcx.append(k)
			convexdecomposedfaces = []
			indices = np.argsort(pcx)
			for x in range(0,len(indices)-1):
				convexdecomposedfaces.append([indices[x],indices[x+1]])
			hull = [[indices[0]],[indices[len(indices)-1]]]					
		else:
			faces = Delaunay(pca_x).simplices
			hull = ConvexHull(pca_x).simplices
			convexdecomposedfaces = iterativeconvexization(faces,len(pca_x[0]),pca_x)
	
		mappedconvexedfaces = []
		for a in convexdecomposedfaces:
			list1 = []
			for y in a:
				if(y>=farthestpoint):
					list1.append(points[y+1])
				else:
					list1.append(points[y])
			mappedconvexedfaces.append(list1)
		for x in hull:
			list1 = []
			list2 = []
			for j in x:
				if(j>=farthestpoint):
					list1.append(points[j+1])
					list2.append(j+1)
				else:
					list1.append(points[j])
					list2.append(j)	
			list1.append(points[farthestpoint])
			list2.append(farthestpoint)
			convexdecomposedfaces.append(list2)
			mappedconvexedfaces.append(list1)
		dimension = len(pca_x[0])
		return polytope(poly,polyunmapped,weight,polytophalfspaces,polyhalfspace,chebR,chebXc,farthestpoint,steriographic_projection,faces,convexdecomposedfaces,mappedconvexedfaces,dimension,pca_x)
	
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
					#print(poly.dimension)
					#print(len(poly.pca_x[0]))
					#print(len(polytop))
					#print(len(poly.pca_x[0])+1)
					if(len(polytop)==poly.dimension):
						mappedconvexfaces = generatesimplexfacets(polytop,poly.dimension)
						convexfaces = generatesimplexfacets([i for i in range(0,len(polytop))],poly.dimension)
						newchild = polytope(polytop,poly.convexdecomposed[k],0,[],[],[],[],-1,[],[],convexfaces,mappedconvexfaces,poly.dimension-1,[])
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
			vertices.append(i)
		polytopaltree.append(vertices)
		dimension1 =dim
		for x in polytopaltree:
			print("polytopes of dimension :: ",dimension1," ::",len(x))
			dimension1 = dimension1-1
	    
		return poly
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
head = tree()  # initialize a polytopla tree
dim= 5  #dimension of a data
datasize = 50
inputpoints = tadasets.dsphere(n=datasize, d=dim-1, r=1, noise=0.1)  # generate d-sphere datatset
points = [i for i in range(0,len(inputpoints))]
computedistancematrix(inputpoints)	#Populate distance matrix
root = head.createheadpolytope(points)  # if data lies in n-dimensions; this will make a surface of a lifted (n+1)-dimensional parbolied its projection will be delaunay triangulation
head.createbootomup(root)  #will need to impliment this function to gro polytopal tree bottom up

#head.display(root)
