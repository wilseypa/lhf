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

def neighbours(polytop,simplices):    #report all the neighbouring polytopes that shares a face with polytope
	polytop = list(polytop)
	adjacency = []
	for y in simplices:
		y = list(y)
		if(polytop != y):
			if(len(intersection(polytop,y))==dim):
				adjacency.append(y)
	return adjacency

def mergeneighbors(polytop,simplices,points):    #merge neigbours simplices that results in maximum convex polytop with neighbors
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


def optimalconvexization(simplices,points):   #Repeate the Maximum convexization for each simplex and returned the sorted list of convex polytopes by weight and vertices
	convexparts = []
	maxdist = []
	averagedist = []
	sizecd = []
	for i in range(0,len(simplices)):
		x= list(simplices[i])
		distance =0
		maxval = 0
		while True:
			x1 = list(mergeneighbors(x,simplices,points))
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

def iterativeconvexization(simplices,dim,points):   #Keep convex decomposition that minimizes convex polytopes required and minimizes maximum weight edge
	valid = []
	convexpartsunion = []
	simplices = list(simplices)
	remainingsimplices = simplices
	premaining = 0
	while(True):
		df = optimalconvexization(remainingsimplices,points)
		i=0
		pointsaddressed = []
		for x in df['convexpart'].tolist():
			check = 1
			for t in convexpartsunion:
				if(len(intersection(t,x)) > dim):
					check = 0
			if(check==1):
				valid.append(i)
				hull = ConvexHull([points[i] for i in x])
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



class polytope:
	def __init__(self,polytope1,weight1,polytophalyspaces,polyhalfspace,chebR1,chebXc1,farthestpoint1,steriographic_projection1,faces,convexdecomposedfaces,dimension,pca_x):
		self.polytopevertices = polytope1  # list of coordinates of polytope vertices
		self.weight = weight1  #weight of polytope
		self.polytophalyspaces = polytophalyspaces  #Hyperplane System of Equations
		self.polyhalfspace = polyhalfspace  # Create polytope object
		self.chebR = chebR1  #polytope incenter
		self.chebXc = chebXc1   #polytope inradius
		self.farthestpoint = farthestpoint1 #polytope fathest point from incenter coordinates
		self.steriographic_projection = steriographic_projection1 # list of corrdinates one dimensional less (Conformal Projection)
		self.faces = faces #list of polytopal faces
		self.convexdecomposedfaces =  convexdecomposedfaces  #PCA Dimensionality Reduction for lower order faces
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
		farthestpoint = []
		steriographic_projection = inputpoints
		faces =  Delaunay(inputpoints).simplices
		convexdecomposedfaces = iterativeconvexization(faces,len(inputpoints[0]),inputpoints)
		dimension = len(inputpoints[0])+1
		pca_x = inputpoints
		return polytope(poly,weight,polytophalyspaces,polyhalfspace,chebR,chebXc,farthestpoint,steriographic_projection,faces,convexdecomposedfaces,dimension,pca_x)

	def createpolytope(self,points,pca_x1):
		poly = points
		print(points)
		polytopal_points = [pca_x1[i] for i in points]
		weight = 0 #delauney retension
		polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopal_points)
		polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
		chebR = polyhalfspace.chebR
		chebXc = polyhalfspace.chebXc
		farthestpoint = fathestpointpolytopeindex(points,chebXc,polytopal_points)
		nearesttofarthest = nearesttofarthestpoint(points,farthestpoint,polytopal_points)
		steriographic_projection = computesteriographicprojection(points ,chebXc,farthestpoint,nearesttofarthest,polytopal_points)
		pca_x = PCAprojection(steriographic_projection)
		faces = Delaunay(pca_x).simplices
		hull = ConvexHull(pca_x).simplices
		convexdecomposedfaces = iterativeconvexization(faces,len(pca_x[0]),pca_x)
		mappedconvexedfaces = []
		for a in convexdecomposedfaces:
			list1 = []
			for y in a:
				list1.append(points[y])
			mappedconvexedfaces.append(list1)
		for x in hull:
			list1 = []
			for j in x:
				list1.append(points[j])
			list1.append(points[farthestpoint])
			mappedconvexedfaces.append(list1)
		print(len(mappedconvexedfaces))
		print(mappedconvexedfaces)
		dimension = len(pca_x[0])+1
		return polytope(poly,weight,polytophalfspaces,polyhalfspace,chebR,chebXc,farthestpoint,steriographic_projection,faces,mappedconvexedfaces,dimension,pca_x)
	
	def createbootomup(self,polytope):  #needs to be implemented
		i = 0
		for polytop in polytope.convexdecomposedfaces:
			newchild = self.createpolytope(polytop,polytope.pca_x)
			print(i)
			i=i+1
			polytope.childrens.append(newchild)
		return polytope
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
dim=3    #dimension of a data
inputpoints = tadasets.dsphere(n=100, d=dim-1, r=1, noise=0.3)  # generate d-sphere datatset
points = [i for i in range(0,len(inputpoints))]
computedistancematrix(inputpoints)	#Populate distance matrix
root = head.createheadpolytope(points)  # if data lies in n-dimensions; this will make a surface of a lifted (n+1)-dimensional parbolied its projection will be delaunay triangulation
head.createbootomup(root)  #will need to impliment this function to gro polytopal tree bottom up

#head.display(root)
