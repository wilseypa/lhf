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
import csv
import fiblat
import itertools
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import math
import sys

maxepsilon = 2
mindist = 0.0001
distancematrix = []	  #initializae distance matrix
delaunayincidence = []
	

def computedistancematrix(points):    #compute Distance Matrix
	global distancematrix    #initializae distance matrix
	distancematrix = []  
	global delaunayincidence
	delaunayincidence = []
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

def neighbours(polytop,simplices,dim1,delaunaypart):    #report all the neighbouring polytopes that shares a face with polytope
	polytop = list(polytop)
	adjacency = []
	for y in simplices:
		y = list(y)
		if(polytop != y):
			if(len(intersection(polytop,y))==dim1):
				adjacency.append(y)
			if(len(intersection(polytop,y))==dim1+1):
				if(y not in delaunaypart):
					delaunaypart.append(y)
	return adjacency,delaunaypart
	
def simplexweight(simplex):
	if(len(simplex)==1):
		return 0
	facets = generatesimplexfacets(simplex,2)
	dist = 0
	for edge in facets:
		if(dist < distanceindex(edge[0],edge[1])):
			dist = distanceindex(edge[0],edge[1])
	return dist

	
def mergeneighbors(polytop,simplices,pointscoord,delaunaypart,points):    #merge neigbours simplices that results in maximum convex polytop with neighbors
	neighbors,delaunaypart = neighbours(polytop,simplices,len(pointscoord[0]),delaunaypart)
	simpl = []
	maxweightedge = 0
	for x in neighbors:
		y = union(x,polytop)
		hull = ConvexHull([pointscoord[points.index(i)] for i in y],qhull_options="QJ")
		hullboundary = {}
		for sublist in hull.simplices:
			for item in sublist:
				if item in hullboundary.keys():
					hullboundary[item] += 1
				else:
					hullboundary[item] = 1
		if(len(list(hullboundary.keys())) == len(y)):
			polytop = union(polytop,x)
			if(maxweightedge<simplexweight(x)):
				maxweightedge = simplexweight(x)
			delaunaypart.append(x)
			
	return list(polytop),delaunaypart,maxweightedge

def getnextsimplex(simplices,x,addresspoints):
	
	addresspoints = union(addresspoints,x)
	for y in simplices:
		if(len(intersection(y,addresspoints)) != len(y)):
			break
	return list(y),addresspoints

def optimalconvexization(simplices,pointscoord,points):   #Repeate the Maximum convexization for each simplex and returned the sorted list of convex polytopes by weight and vertices
	convexparts = []
	maxdist = []
	points= sorted(points)
	averagedist = []
	totalpoints = []
	[totalpoints.append(k) for y in simplices for k in y if k not in totalpoints]
	delaunayparts = []
	sizecd = []
	updatedsimplices = simplices
	addressedpoints = []
	x = []
	for i in range(0,len(simplices)):
		x,addressedpoints = getnextsimplex(simplices,x,addressedpoints)
		delaunaypart = []
		delaunaypart.append(x)
		distance =0
		maxweight = simplexweight(x)
		maxval = 0
		while True:
			x1,delaunaypart,w = mergeneighbors(x,simplices,pointscoord,delaunaypart,points)
			if(maxweight<w):
				maxweight = w
			if(x1==x):
				break
			x = x1
		if x not in convexparts:
			maxval=0
			distance = 0
			convexparts.append(x)
			maxdist.append(maxweight)
			delaunayparts.append(delaunaypart)
			hull = ConvexHull([pointscoord[points.index(i)] for i in x],qhull_options="QJ")
			for p in hull.simplices:
				distance = distance + math.dist(pointscoord[points.index(x[p[0]])],pointscoord[points.index(x[p[1]])])
				if(maxval < math.dist(pointscoord[points.index(x[p[0]])],pointscoord[points.index(x[p[1]])])):
					maxval = math.dist(pointscoord[points.index(x[p[0]])],pointscoord[points.index(x[p[1]])])
			maxdist.append(maxval)
			averagedist.append(distance)
			sizecd.append(len(x))
		k=set()
		for tp in convexparts:
			for lm in tp:
				k.add(lm)
		if(len(k)==len(totalpoints)):
			break;
	df = pd.DataFrame(list(zip(convexparts,delaunayparts,maxdist,averagedist,sizecd,maxdist)),columns =['convexpart','delaunayparts','maxdist','averagedist','sizecd','weight'])
	df = df.sort_values(by = 'sizecd',ascending = False)
	df = df.sort_values(by = 'maxdist',ascending = True)
	return df

def iterativeconvexization(simplices,dim,pointscoord):   #Keep convex decomposition that minimizes convex polytopes required and minimizes maximum weight edge
	valid = []
	convexpartsunion = []
	delaunaypartsfinal = []
	weights = []
	points = []
	[points.append(i) for k in simplices for i in k if i not in points]
	delaunayparts = []
	simplices = list(simplices)
	remainingsimplices = simplices
	premaining = 0
	while(True):
		df = optimalconvexization(remainingsimplices,pointscoord,points)
		i=0
		pointsaddressed = []
		for x,y,weight in zip(df['convexpart'].tolist(),df['delaunayparts'].tolist(),df['weight'].tolist()):
			check = 1
			for t in convexpartsunion:
				if(len(intersection(t,x)) > dim):
					check = 0
			if(check==1):
				valid.append(i)
				hull = ConvexHull([pointscoord[points.index(i)] for i in x],qhull_options="QJ")
				pointsaddressed = union(pointsaddressed,x)
				convexpartsunion.append(x)
				delaunayparts.append(y)
				weights.append(weight)
				delaunaypartsfinal.append(y)
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
					weights.append(simplexweight(r))
					delaunayparts.append([list(r)])
					delaunaypartsfinal.append([list(r)])
			break
		premaining = remaining
	convexpartssorted = []
	for x in convexpartsunion:
		convexpartssorted.append(sorted(x))
	return convexpartssorted,delaunayparts,weights


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

def mapdelaunayindices(fromindices,toindices):
	newindices = []
	for x in fromindices:
		lst = []
		for y in x:
			l = []
			for z in y:
				l.append(toindices[z])
			lst.append(l)
		newindices.append(lst)
	return np.array(newindices,dtype=object)
	

def hyperplane(points):
   X=np.matrix(points)
   k=np.ones((len(points[0]),1))
   a=np.matrix.dot(np.linalg.pinv(X), k)
   
   return a
 
def solutionlinalg(matrixA,matrixB):
	return np.linalg.pinv(matrixA)*matrixB  

def PCAprojection(pts):
	pca = PCA(n_components=len(pts[0])-1).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x
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
	return (u - proj_of_u_on_n)+planecentroid
'''	
def maximumfacetprojection(coordinates_points):
	polytopalpoints = coordinates_points
	
	polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopalpoints)
	polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
	chebR = polyhalfspace.chebR
	chebXc = polyhalfspace.chebXc
	projectedpts = findfarthestfacetprojection(coordinates_points)
	triangulation = ConvexHull(coordinates_points).simplices
	convexdecomposition = iterativeconvexization(triangulation,len(projectedpts[0]),projectedpts)
	maxleng = 0
	for x in convexdecomposition:
		if(maxleng<len(x)):
			maxfacet = x
			maxleng = len(x)
	ctr, normal  = planeFit([coordinates_points[k] for k in maxfacet])
	projectedpoints = []
	for x in maxfacet:
		projectedpoints.append(projectonplane(coordinates_points[x],normal,ctr))
	print(maxfacet)	
	newcoordinates = []
	k=0
	newhull = ConvexHull(coordinates_points).simplices
	print(len(coordinates_points)," ",len(coordinates_points))
	print(newhull)
	updatedhull = []
	updatedhull.append(np.array(maxfacet))
	for tri in newhull:
		if(all(x in maxfacet for x in tri)):
			continue
		else:
			updatedhull.append(tri)
	updatedhull = np.array(updatedhull)
	print(updatedhull)
	mindist = 100000
	for x in maxfacet:
		dist  = math.dist(chebXc,coordinates_points[x])
		if(dist<mindist):
			mindist = dist
	for x in range(0,len(maxfacet)):
		dist  = math.dist(chebXc,projectedpoints[x])
		if(dist<mindist):
			mindist = dist

	multiplier = mindist/math.dist(chebXc,ctr)
	stereoprojection_point = []
	for x,y in zip(ctr,chebXc):
		stereoprojection_point.append(multiplier*(x-y)+y)
	dim = len(coordinates_points[0])
	stereoprojection_plane = hyperplane([projectedpoints[t] for t in range(0,dim)])
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
	return PCAprojection(pts),updatedhull
'''	
def findfarthestfacetprojection(coordinates_points,farthestface):  #compute farthest facet steriographic projection
	polytopalpoints = coordinates_points
	polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopalpoints)
	polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
	chebR = polyhalfspace.chebR
	chebXc = polyhalfspace.chebXc
	faces = ConvexHull(polytopalpoints,qhull_options="QJ").simplices
	farthest = 0
	normal = True
	if(len(farthestface) == len(coordinates_points[0])):
		normal = False
		
	if(farthestface == []):
		farthestface = faces[0]
		normal = False
	centroid1 = []
	projectionfacet = []
	if(normal):
		farthestfacet = []
		hullfacet = ConvexHull([coordinates_points[i] for i in farthestface],qhull_options="QJ").simplices
		hullfacetmappeds = mapindices(hullfacet,farthestface)
		hullfacetmapped = []
		[hullfacetmapped.append(sorted(g)) for g in hullfacetmappeds]
		dst = 99999
		for x in hullfacetmapped:
			ponts = [coordinates_points[y] for y in x]
			cntr = sum(np.transpose(list(list(y) for y in list(np.transpose(ponts)))))/len(ponts[0])
			if(math.dist(chebXc,cntr)<dst):
				farthestfacet = x
				dst = math.dist(chebXc,cntr)
		ponts = [coordinates_points[y] for y in farthestface]
		stereoprojection_point = sum(np.transpose(list(list(y) for y in list(np.transpose(ponts)))))/len(ponts[0])
		stereoprojection_plane = hyperplane([coordinates_points[t] for t in farthestfacet])
		cofficient = [y for j in np.array(stereoprojection_plane) for y in np.array(j)]
		constant  = -1
		pts = []
		count = 0
		for i in range(0,len(coordinates_points)):
			if(i not in farthestface):
				count = count+1
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
			else:
				done = False
				for x in faces:
					if(done):
						break
					if(len(intersection(x,[i]))==1):
						if(len(intersection(x,farthestface))==1):
							count = count+1
							facetlower = [coordinates_points[y] for y in x if y not in [i]]
							centroid = sum(np.transpose(list(list(x) for x in list(np.transpose(facetlower)))))/len(facetlower[0])
							diff = centroid - coordinates_points[i]
							constterm =0;
							diffterm = 0;
							for x in range(0,len(cofficient)):
								constterm+= coordinates_points[i][x]*cofficient[x]
								diffterm+= diff[x]*cofficient[x]
							if(diffterm!=0):
								t = ((-1)*constant - constterm)/diffterm
								pnt = t*(centroid-coordinates_points[i])+coordinates_points[i]
								pts.append(pnt)
							done = True
				if(done == False):
					for x in faces:
						if(done):
							break
						notfarthestface = []
						[notfarthestface.append(lg) for lg in range(0,len(coordinates_points)) if lg not in farthestface]
						if(len(intersection(x,[i]))==1 and len(intersection(x,notfarthestface))==1):
							count = count+1
							upperfacet = [coordinates_points[y] for y in x if y not in notfarthestface]
							centroid = sum(np.transpose(list(list(x) for x in list(np.transpose(upperfacet)))))/len(upperfacet[0])
							diff = centroid - coordinates_points[i]
							constterm =0;
							diffterm = 0;
							for x in range(0,len(cofficient)):
								constterm+= coordinates_points[i][x]*cofficient[x]
								diffterm+= diff[x]*cofficient[x]
							if(diffterm!=0):
								t = ((-1)*constant - constterm)/diffterm
								pnt = t*(centroid-coordinates_points[i])+coordinates_points[i]
								pts.append(pnt)
							done = True
				
						
		return PCAprojection(pts)
	else:					
		'''
		for edges in faces:
			pts = [coordinates_points[x] for x in edges]
			centroid = sum(np.transpose(list(list(x) for x in list(np.transpose(pts)))))/len(pts[0])
			if(math.dist(centroid,chebXc)>farthest):
				farthest = math.dist(centroid,chebXc)
				farthestface = edges
				centroid1 = centroid
		farthestface = list(farthestface)
		'''
		pts = [coordinates_points[x] for x in farthestface]
		centroid1 = sum(np.transpose(list(list(x) for x in list(np.transpose(pts)))))/len(pts[0])
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
		constant  = -1
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
def hullfromtriangulation(simplices):
	hull= []
	for x in simplices:
		facets = generatesimplexfacets(x,len(x)-1)
		for f in facets:
			valid = True
			for y in simplices:
				if(valid==False):
					break
				if(x !=y):
					if(len(intersection(f,y))==len(f)):
						valid=False
						break
			if(valid==True):
				if(sorted(f) not in hull):
					hull.append(sorted(f))
	return hull
def generatelowerorderdecompositions(indices,coordinates):
	level = []
	indices = indices
	coords = coordinates
	points = []
	indexs = []
	[points.append(x) for z in indices for x in z if x not in points]
	[indexs.append(y) for y in range(0,len(points))]
	points = sorted(points)
	triangulation = hullfromtriangulation(indices) 
	projectedpts = findfarthestfacetprojection(coords,[])
	convexdecomposition,delaunayparts,weights = iterativeconvexization(triangulation,len(projectedpts[0]),projectedpts)
	convexhullpoints = ConvexHull(projectedpts).simplices
	mappedhullpoints = mapindices(convexhullpoints,points)
	hullpoints = []
	[hullpoints.append(x) for z in mappedhullpoints for x in z if x not in hullpoints]
	currentfacetlen = 0
	projectionfacet = []
	for x in convexdecomposition:
		if(len(intersection(x,hullpoints))==0):
			if(currentfacetlen<len(x)):
				projectionfacet = x
				currentfacetlen = len(x)
	projectionfacetmapped = []
	[projectionfacetmapped.append(points.index(x)) for x in projectionfacet]
	projectedpts = findfarthestfacetprojection(coords,projectionfacetmapped)
	convexdecomposition,delaunayparts,weights = iterativeconvexization(triangulation,len(projectedpts[0]),projectedpts)
	convexdecompositionmapped = convexdecomposition
	delaunaypartsmapped = delaunayparts
	for (x,y,z,w) in zip(convexdecomposition,convexdecompositionmapped,delaunaypartsmapped,weights):
		x = sorted(x)
		lowercoordinates = [projectedpts[points.index(i)] for i in x]
		level.append([z,lowercoordinates,w])
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
			triangulation = ConvexHull(coords,qhull_options="QJ").simplices
			mappedindices = mapindices(triangulation,points)
			[level.append(sorted(t)) for t in mappedindices.tolist() if sorted(t) not in level]
		else:
			facets = generatesimplexfacets(points,2)
			[level.append(sorted(t)) for t in facets if sorted(t) not in level]
	polytopaltree.append(level)
	print("Polytopes at Dimension ::",1,"::",len(level))
	level = []
	for x in range(0,len(inputpoints)):
		level.append([x])
	polytopaltree.append(level)
	print("Polytopes at Dimension ::",0,"::",len(level))
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
	convexdecomposition,delaunayparts = iterativeconvexization(triangulation,len(coordinates[0]),coordinates)
	level = []
	for (x,z) in zip(convexdecomposition,delaunayparts):
		x = sorted(x)
		lowercoordinates = [coordinates[i] for i in x]
		level.append([z,lowercoordinates])
	print("Polytopes at Dimension ::",dim)
	#cascade down to Triangles
	for d in range(1,dim-1):
		convexdecompositions = copy.deepcopy(polytopaltree[d-1])
		level = []
		levelcheck = []
		lenthr = 0
		for cd in convexdecompositions:
			if(len(cd[0])!=1):
				lowercd = ConvexHull(cd[1],qhull_options="QJ").simplices
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
		print("Polytopes at Dimension ::",(dim-d))
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
	convexdecomposition,delaunayparts,weights = iterativeconvexization(triangulation,len(coordinates[0]),coordinates)
	level = []
	for x,z,w in zip(convexdecomposition,delaunayparts,weights):
		x = sorted(x)
		lowercoordinates = [coordinates[i] for i in x]
		level.append([z,lowercoordinates,w])
	polytopaltree.append(level)
	print("Polytopes at Dimension ::",dim,"::",len(level))
	#cascade down to Triangles
	for d in range(1,dim):
		convexdecompositions = copy.deepcopy(polytopaltree[d-1])
		level = []
		levelcheck = []
		lenthr = 0
		for cd in convexdecompositions:
			if(len(cd[0])!=1):
				if(d == dim-1):
					facets = ConvexHull(cd[1],qhull_options="QJ").simplices
					for x in facets:
						y = sorted(x)
						y= [y]
						if y not in levelcheck:
							levelcheck.append(y)
							level.append([y,y,simplexweight(y)])
					lenthr = lenthr + len(facets)
				else:
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
						level.append([y,y,simplexweight(y)])
				lenthr = lenthr + len(facets)
		polytopaltree.append(level)
		print("Polytopes at Dimension ::",(dim-d),"::",len(level))
	level = []
	for x in range(0,len(inputpoints)):
		level.append([[[x]]])
	polytopaltree.append(level)
	print("Polytopes at Dimension ::",0,"::",len(level))
	#Add vertices
	#polytopaltree = addEdgesAndVertices(polytopaltree)
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
			if(self.polytop==nxt.polytop):
				return True
			for (x,y) in zip(reversed(self.polytop),reversed(nxt.polytop)):
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
'''
def earcliptriangulate(inputpoints,coordinates):
	points = []
	[points.append(x) for y in inputpoints for x in y if x not in points]
	hull = ConvexHull(coordinates).simplices
	dim = len(hull[0])
	triangulation = []
	polytopalpoints = coordinates
	
	polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopalpoints)
	polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
	chebR = polyhalfspace.chebR
	chebXc = polyhalfspace.chebXc
	first = True
	removed = []
	while(len(polytopalpoints)>dim+1):
		polytophalfspaces = polyfun.compute_polytope_halfspaces(polytopalpoints)
		polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
		chebR = polyhalfspace.chebR
		chebXc = polyhalfspace.chebXc
		farthestpoint = 0
		maxdist = 0
		faces = ConvexHull(polytopalpoints).simplices
		print(faces)
		if(first==False):
			print(faces)
			print(cntpoints)
			for fc in faces:
				if(all(x in cntpoints for x in fc)):
					lst = list(fc)
					print(fc)
					lst.append(farthestpoint)
					triangulation.append(lst)
					print(triangulation)
					input()
		for x in range(0,len(polytopalpoints)):
			distance = math.dist(chebXc,polytopalpoints[x])
			if(maxdist<distance):
				maxdist = distance
				farthestpoint = x
		connectedpoints = []
		for fc in faces:
			if farthestpoint in fc:
				connectedpoints.append(fc)
		cntpoints1 = []
		[cntpoints1.append(x) for y in connectedpoints for x in y if x not in cntpoints1 if x !=farthestpoint]
		cntpoints = []
		for x in cntpoints1:
			if(x>farthestpoint):
				cntpoints.append(x-1)
			else:
				cntpoints.append(x)
		polytopalpoints.pop(farthestpoint) 
		removed.append(farthestpoint)
		removed = sorted(removed)
		first = False
		print("here")
		print(cntpoints)
		print(faces)
		print(farthestpoint)
		print(connectedpoints)
	#faces = ConvexHull(polytopalpoints).simplices
	#farthest = 0
	#farthestface = []
	input()
	return triangulation;
'''
def assignweightsold(polytopaltree):
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
							#if(incidenceindexreturn(edge[0],edge[1])==1):
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

def assignweights(polytopaltree):
	orderedpolytoparraylist = []
	for dimensionlist in polytopaltree:
		dimensionwiseorderedpolytopes = set()
		for polytope in dimensionlist:	
			maxedge = 0
			polytp = []
			for poly in polytope[0]:
				if(maxedge<simplexweight(poly)):
					maxedge = simplexweight(poly)
				polytp.append([poly,simplexweight(poly)])
			node = orderedarraylistnode(polytp,maxedge)
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
			bettieTable.append(tuple(bettitableentry))
			if(mstSize >= len(inputpoints)-1):
				for i in range(0,len(inputpoints)):
					if(ds.find(i) == i):
						bettitableentry = [0,0,maxepsilon]
						bettieTable.append(tuple(bettitableentry))
				return pivots
	return pivots

def getDimEdges(dimension):
	if(dimension==0 or dimension == dim):
		edg =  polytopaltree[dim -dimension]
	else:
		edg = polytopaltree[dim -dimension-1]
		dimensionwiseorderedpolytopes = set()
		for x in edg:
			boundary = generatesimplexfacets(x.polytop,len(x.polytop)-1)
			for poly in boundary:
				localmax =0
				edges = generatesimplexfacets(poly,2)
				for edge in edges:
					if(localmax<distanceindex(edge[0],edge[1])):
							localmax = distanceindex(edge[0],edge[1])
				node = orderedarraylistnodeorig(poly,localmax)
				dimensionwiseorderedpolytopes.add(node)
		sorted_list = list(dimensionwiseorderedpolytopes)
		sorted_list.sort(key = letter_cmp_key)
		edg = sorted_list
	return edg
	
def getAllCofacets(e,dimension):    #will update this function a little
	returnlist = []
	for x in polytopaltree[dim-dimension-1]:
		if(len(intersection(x.polytop,e.polytop)) == len(e.polytop)):
			returnlist.append(x)
	return returnlist

def persistenceByDimension( edges, pivots, dimension):
	pivotcounter = len(pivots)-1
	edges.sort(key = letter_cmp_key)
	nextPivots = []	
	'''
	print("Pivots")
	for x in pivots:
		print(x.polytop,"  ",x.weight)
	print("simpklices")
	for y in edges:
		print(y.polytop," ",y.weight)
	pivots.sort(key = letter_cmp_key)
	input()
	'''
	v = {} 
	pivotPairs = {}
	for e in edges:
		if(pivotcounter < 0 or pivots[pivotcounter].polytop != e.polytop):
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
					'''
					print("Facelist1")
					for x in faceList:
						print(x.polytop,x.weight)
					print("Facelist2")
					'''
					pivot =	hq.heappop(faceList)
					'''
					print(pivot.polytop,pivot.weight)
					if(faceList != []):
						print(faceList[0].polytop,faceList[0].weight)
					'''				
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

def fastpersistance(polytopalcomplex):
	vertices = getDimEdges(0)
	edges = getDimEdges(1)
	pivots  = minimumspanningtree(edges)
	
	for d  in range(1,dim):
		if(d != 1):
			 edges = getDimEdges(d)
		pivots = persistenceByDimension(edges, pivots, d)
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
def disintegrateold(polytopaltree):
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
	
def disintegrate(polytopaltree):
	orderedpolytoparraylist = []
	for dimensionlist in polytopaltree:
		dimensionwiseorderedpolytopes = set()
		for polytope in dimensionlist:
			for part in polytope.polytopparts:
				node = orderedarraylistnodeorig(part[0],part[1])
				dimensionwiseorderedpolytopes.add(node)
		sorted_list = list(dimensionwiseorderedpolytopes)
		sorted_list.sort(key = letter_cmp_key)
		orderedpolytoparraylist.append(sorted_list)
	return orderedpolytoparraylist



def checkplotsimplices(simplices):
	newpolyparts = []
	indices = []
	[indices.append(i) for x in simplices for i in x if i not in indices]
	for x in simplices:
		for t in ConvexHull([inputpoints[i] for i in x],qhull_options="QJ").simplices:
			t = [x[i] for i in t]
			newpolyparts.append(t)
	polyparts = newpolyparts
	hull = polyparts#ConvexHull(points).simplices
	newlist = []
	i = 1
	fig = plt.figure()
	ax = Axes3D(fig, auto_add_to_figure=False)
	fig.add_axes(ax)	
	for y in hull:
		#newlist.append(y)
		for x in [y]:
			pts = [inputpoints[x1] for x1 in x]
			ptrans = np.transpose(pts)
			verts = [list(zip(ptrans[0],ptrans[1],ptrans[2]))]
			ax.scatter3D(ptrans[0],ptrans[1],ptrans[2])
			ax.add_collection3d(Poly3DCollection(verts,facecolors='b', edgecolor = 'black', linewidths=1, alpha=0.15))
		if(i%12==0):
			plt.show()
		i=i+1

def plotpolytop(polyparts,d):
	if(d==3):
		newpolyparts = []
		for x in polyparts:
			x=x[0]
			for t in ConvexHull([inputpoints[i] for i in x],qhull_options="QJ").simplices:
				t = [x[i] for i in t]
				if(sorted(t) not in newpolyparts):
					newpolyparts.append(sorted(t))
		polyparts = newpolyparts
	hull = polyparts#ConvexHull(points).simplices
	fig = plt.figure()
	ax = Axes3D(fig, auto_add_to_figure=False)
	fig.add_axes(ax)
	for x in hull:
		if(d!=3):
			#x=x[0]
			x=x
		pts = [inputpoints[x1] for x1 in x]
		ptrans = np.transpose(pts)
		verts = [list(zip(ptrans[0],ptrans[1],ptrans[2]))]
		ax.scatter3D(ptrans[0],ptrans[1],ptrans[2])
		ax.add_collection3d(Poly3DCollection(verts,facecolors='b', edgecolor = 'black', linewidths=1, alpha=0.15))
	plt.show()

def genPermutahedron(a, size):
	global ptspermutahedron
	if size == 1:
		ptspermutahedron.append(np.array(a))
		return
	for i in range(size):
		genPermutahedron(a, size-1)
		if size & 1:
			a[0], a[size-1] = a[size-1], a[0]
		else:
			a[i], a[size-1] = a[size-1], a[i]
def swapdifferbyone(pnt1,pnt2):
	cnt=0
	lst = []
	for i in range(0,len(pnt1)):
		if(pnt1[i]!=pnt2[i]):
			cnt=cnt+1
			lst.append(pnt1[i])
	if(cnt==2 and abs(lst[0]-lst[1])==1):
		return True
	else:
		return False
		
def permutahedronincidence(points):
	incidence = []
	for i in range(0,len(points)):
		lst = []
		for j in range(0,len(points)):
			if(swapdifferbyone(points[i],points[j])):
				lst.append(1)
			else:
				lst.append(0)
		incidence.append(lst)
	return incidence

def generatehypercube(dimensions):
	result = list(itertools.product('10', repeat=dimensions))
	inputpoints = [[int(y) for y in x] for x in result]
	return np.array(inputpoints)-0.5

def differbyone(pnt1,pnt2):
	cnt=0
	for i in range(0,len(pnt1)):
		if(pnt1[i]!=pnt2[i]):
			cnt=cnt+1
	if(cnt>1):
		return False
	else:
		return True

def hypercubeincidence(points):
	incidence = []
	for i in range(0,len(points)):
		lst = []
		for j in range(0,len(points)):
			if(differbyone(points[i],points[j])):
				lst.append(1)
			else:
				lst.append(0)
		incidence.append(lst)
	return incidence
def refine(points):
	rmlist = []
	for i in points:
		for j in points:
			if((i!=j).all()):
				if(math.dist(i,j)<mindist):
					rmlist.append(j)
	for i in rmlist:
		points.pop(i)
	return points

def PCAprojection(pts):
	pca = PCA(n_components=len(pts[0])-1).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x

def generatesimplexfacets(poly,dimen):  #generate simplices of dimension dimen
	tp = list(combinations(poly, dimen))
	listreturn = []
	for x in tp:
		listreturn.append(list(x))
	return listreturn

def addgaussiannoise(unitsphere,noiseper):
	x, y = np.meshgrid(np.linspace(-1,1,len(unitsphere[0])), np.linspace(-1,1,len(unitsphere)))
	dst = np.sqrt(x*x+y*y)
	sigma = 1
	muu = 0.0001
 	# Calculating Gaussian array
	gauss = np.exp(-( (dst-muu)**2 / ( 2.0 * sigma**2 ) ) )
	unitsphere = np.array(unitsphere) + gauss*noiseper
	return unitsphere
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

def getstarted(datatype,dimension,datasize,noiseper):
	unitsphere = []
	filenameinput = 'input'+datatype
	filenameoutput = 'output'+datatype
	if(datatype=='0'):
		print("Please Choose Data Type::")
		print("1: Hypertetrahedron Sphere")
		print("2: HypercubeSphere Sphere")
		print("3: Hyperpermutahedron Sphere")
		print("4: Hyperfibonnacci Sphere")
		print("5: Read From File")
		datatype = input()
	if(datatype!='5'):
		if(dimension=='-1'):
			print("Dimension ::")
			dimension = input()
		filenameinput = filenameinput+'dim'+dimension
		filenameoutput = filenameoutput+'dim'+dimension
		dim = int(dimension)
		if(datasize=='-1'):
			print("Intended DataSize ::")
			datasize = input()
		filenameinput = filenameinput+'size'+datasize
		filenameoutput = filenameoutput+'size'+datasize
		datasize = int(datasize)
		if(noiseper=='-1'):
			print("Gaussian noise in percentage ::")
			noiseper = input()
		filenameinput = filenameinput+'noise'+ noiseper
		filenameoutput = filenameoutput+'noise'+ noiseper
		noiseper = int(noiseper)/1000
	if datatype=='1':  #Tetrahedron
		edgecount = (int)(dim*(dim+1)/2)
		points = max(datasize,dim+1)
		pointperedge = (int)(points/edgecount) +1
		poly = [i for i in range(0,dim+1)]	
        
		for dimension in range(dim,dim+1):
			points = []
			pnts = dimension+1
			cnt  =0
			origin = [0 for x in range(0,dimension)]
			for x in range(0,pnts):
				coords = []
				for y in range(0,dimension+1):
					if(y==cnt):
						coords.append(1)
					else:
						coords.append(0)
				points.append(coords)
				cnt = cnt+1
			newpoints = PCAprojection(points)
			points = set()
			edges = generatesimplexfacets(poly,2)
			for e in edges:
				for t in np.linspace(0.1, 0.9, pointperedge):
					points.add(tuple(t*(newpoints[e[0]])-newpoints[e[1]]*(t-1)))
			for x in points:
				x = np.array(x)
				unitsphere.append(x/math.dist(origin,x))
	elif datatype=='2':	#"HyperCube"
		vertices = 2**dim
		edgecount = 2**(dim-1)*dim
		points = max(datasize,vertices)
		pointperedge = (int)(points/edgecount) +1
		poly = [i for i in range(0,vertices)]
		incidence = []
		for dimension in range(dim,dim+1):
			pnts = dimension+1
			cnt  =0
			origin = [0 for x in range(0,dimension)]
			newpoints = generatehypercube(dim)
			incidence = hypercubeincidence(newpoints)
			points = set()
			edges = generatesimplexfacets(poly,2)
			for e in edges:
				if(incidence[e[0]][e[1]]==1):
					for t in np.linspace(0.1, 0.9, pointperedge):
						points.add(tuple(t*(newpoints[e[0]])-newpoints[e[1]]*(t-1)))
			for x in points:
				x = np.array(x)
				unitsphere.append(x/math.dist(origin,x))
	elif datatype=='3':	#Permutahedron
		dim = dim+1
		b = [i for i in range(0,dim)]
		vertices = math.factorial(dim)
		edgecount = ((dim-1)*math.factorial(dim))/2
		points = max(datasize,vertices)
		pointperedge = (int)(points/edgecount) +1
		poly = [i for i in range(0,vertices)]
		incidence = []
		for dimension in range(dim-1,dim):
			pnts = dimension+1
			cnt  =0
			origin = [0 for x in range(0,dimension)]
			size = len(b)
			global ptspermutahedron
			ptspermutahedron = []
			genPermutahedron(b,size)
			incidence = permutahedronincidence(ptspermutahedron)
			pca_x = PCAprojection(ptspermutahedron)
			pts = pca_x
			pts = np.array(pts)
			newpoints = pts
			points = set()
			edges = generatesimplexfacets(poly,2)
			for e in edges:
				if(incidence[e[0]][e[1]]==1):
					for t in np.linspace(0.1, 0.9, pointperedge):
						points.add(tuple(t*(newpoints[e[0]])-newpoints[e[1]]*(t-1)))
			for x in points:
				x = np.array(x)
				unitsphere.append(x/math.dist(origin,x))
				
	elif datatype=='4': #FibbonacciLaticce
		unitsphere =  np.array(fiblat.sphere_lattice(dim,datasize))
	
	elif datatype=='5': #ReadFile
		print("Enter File Name")
		filename = input()
		with open(filaname, "r") as f:
			reader = csv.reader(f)
			for line in reader:
				k=[]
				for x in line:
					k.append(float(x))
				unitsphere.append(k)
	
	if(datatype!='5'):
		unitsphere  = refine(unitsphere)
		unitsphere = addgaussiannoise(unitsphere,noiseper)
		print("Final Data Size::",datasize," in dimension ::",dim)
		with open(filenameinput, "w", newline="") as f:
			writer = csv.writer(f)
			writer.writerows(unitsphere)
	return unitsphere,filenameoutput
ptspermutahedron = []
datatype = [1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4]
dimension = [3,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5,3,3,3,3,4,4,4,4,5,5,5,5,3,3,3,3,4,4,4,4,5,5,5,5,3,3,3,3,4,4,4,4,5,5,5,5,3,3,3,3,4,4,4,4,5,5,5,5]
datasize = [100,100,100,100,250,250,250,250,300,300,300,300,400,400,400,400,200,200,200,200,300,300,300,300,400,400,400,400,200,200,200,200,300,300,300,300,400,400,400,400,200,200,200,200,300,300,300,300,400,400,400,400,200,200,200,200,300,300,300,300,400,400,400,400]
noiseper = [0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100]
for d,dim,s,n in zip(datatype,dimension,datasize,noiseper):
	print(" Type:: ",d, " Dimension :: ",dim," Size :: ",s," Noise:: ",n)
	inputpoints,outputfilename = getstarted(str(d),str(dim),str(s),str(n))
	datasize = len(inputpoints)
	dim=len(inputpoints[0])
	'''
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	zdata = [item[0] for item in inputpoints]
	xdata = [item[1] for item in inputpoints]
	ydata = [item[2] for item in inputpoints]
	ax.scatter3D(xdata, ydata, zdata)
	plt.show()
	'''
	letter_cmp_key = cmp_to_key(letter_cmp)    #comparison function
	polytopaltree = createpolytopaltree(inputpoints)
	polytopaltree = assignweights(polytopaltree)
	polytopaltree = disintegrate(polytopaltree)
	
	y = dim
	for x in polytopaltree:
		print("Simplices of Dimension :: ",y,"::--------",len(x))
		y=y-1

	bettieTable = []
	fastpersistance(polytopaltree)
	dimcount=[0 for i in range(0,dim)]
	table = [[] for i in range(0,dim)]
	for x in bettieTable:
		table[x[0]].append(x)
		dimcount[x[0]]= dimcount[x[0]]+1
	for i in range(0,dim):
		print("Dimemnsion ",i," Betti Count ::",dimcount[i])
	writebetties(table,outputfilename)




'''
dimcount=[0 for i in range(0,dim)]

table = [[] for i in range(0,dim)]

for x in bettieTable:
	table[x[0]].append(x)
	dimcount[x[0]]= dimcount[x[0]]+1
for i in range(0,dim):
	print("Dimemnsion ",i," Betti Count ::",dimcount[i])

colors = ["Orange","yellow","Green",'Blue','Red','Black',"pink"]

for d in range(0,dim):
	birth = [x[1] for x in table[d]]
	death = [x[2] for x in table[d]]
	df = pd.DataFrame(list(zip(birth,death)),columns =['birth','death'])
	df = df.sort_values(by = 'death',ascending = True)
	plt.scatter(df["birth"], df["death"],color=colors[d],s=6)
plt.axline([0, 0], [2, 2],linewidth=1,color="black")
plt.savefig("outputBCPolytopal.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()

for d in range(0,dim):
	birth = [x[1] for x in table[d]]
	death = [x[2] for x in table[d]]
	df = pd.DataFrame(list(zip(birth,death)),columns =['birth','death'])
	df = df.sort_values(by = 'death',ascending = True)
	counter = []
	[counter.append(j+i) for j in range(0,len(table[d]))]
	i=i+len(table[d])
	plt.plot([df["birth"], df["death"]], [counter, counter],color=colors[d],linestyle='solid',linewidth=1)
plt.show()
plt.savefig("outputPIpolytopal.pdf", bbox_inches = 'tight',pad_inches = 0)
'''
