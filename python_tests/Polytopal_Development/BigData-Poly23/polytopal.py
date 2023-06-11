import tadasets
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from itertools import combinations
import pypoman.duality as polyfun
import polytope as pc
import math
import fiblat
from sklearn.decomposition import PCA
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from scipy.spatial.transform import Rotation as R
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from numpy.linalg import norm
import sys
sys.path.insert(0, '/home/rohit/Desktop/lhfModified/lhf/python_tests/Polytopal_Development')
from FinalPolytopalImplementation import iterativeconvexization,optimalconvexization,mergeneighbors,neighbours,incidenceindexreturn,setincidenceindex,intersection,union,computedistancematrix,plotpolytop
 
def define_circle(p1, p2, p3):
    """
    Returns the center and radius of the circle passing the given 3 points.
    In case the 3 points form a line, returns (None, infinity).
    """
    temp = p2[0] * p2[0] + p2[1] * p2[1]
    bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
    cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
    det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])
    
    if abs(det) < 1.0e-6:
        return (None, np.inf)
    
    # Center of circle
    cx = (bc*(p2[1] - p3[1]) - cd*(p1[1] - p2[1])) / det
    cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det
    
    radius = np.sqrt((cx - p1[0])**2 + (cy - p1[1])**2)
    return ((cx, cy), radius)
    

def rotation(a,b):
	"""
	rotate vector b to a
	"""
	a,b = a/ np.linalg.norm(a), b/ np.linalg.norm(b)
	v = np.zeros((len(a),len(a)))
	for i in range(len(a)):
		for j in range(len(b)):
			v[i,j] = a[i]*b[j]-a[j]*b[i]
	c = np.dot(a,b)
	I = np.identity(len(a))
	if(c==-1):
		r = I + v
	else:
		r = I + v + np.matmul(v,v) * (1/(1+c))
	return r

def simplexDiagonalWithCenter(radius,dim):
	"""
	Generate Diagonal of simplex
	"""
	diags = np.zeros((dim,dim))
	for i in range(dim):
		diags[i][i] = 1
	pts = []
	for x in diags:
		x = np.array(x)
		pts.append((x*radius)/np.linalg.norm(x))
	pts.append([0 for i in range(dim)])
	pts = np.array(pts)
	points = PCAprojection(pts)
	data = np.c_[points, np.zeros(len(points))]
	return data


def simplexDiagonals(radius,dim):
	"""
	Generate Diagonal of simplex in lower dimension
	"""
	diags = np.zeros((dim,dim))
	for i in range(dim):
		diags[i][i] = 1
	pts = []
	for x in diags:
		x = np.array(x)
		pts.append((x*radius)/np.linalg.norm(x))
	pts = np.array(pts)
	points = PCAprojection(pts)
	return points

	
def intersection(lst1, lst2):  
	"""
	return intersection of two lists
	"""
	return list(set(lst1) & set(lst2))

def PCAprojection(pts):
	"""
	PCA Projection
	"""
	pca = PCA(n_components=len(pts[0])-1).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x

def PCAprojectiondim(pts,dim):
	"""
	PCA projection also Returns the pcs object
	"""
	pca = PCA(n_components=dim).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x,pca


def projectonplane(p1,p2,points1):
	"""
	project points on plane orthogonal to (p1-p2) and having point p1
	"""
	cofficient = p1 - p2
	constant  = 0
	point_on_plane = p1
	for a,b in zip(point_on_plane,cofficient):
		constant = constant+a*b
	constant = (-1)*constant 
	'''
	x = np.linspace(-100, 100, 100)
	y = np.linspace(-100, 100, 100)

	x, y = np.meshgrid(x, y)
	eq = -1*(cofficient[0] * x + cofficient[1] * y + constant)/cofficient[2]
	ax.plot_surface(x, y, eq,color="cyan")
	'''
	pts1 = []
	for pt in points1:
		diff = pt-p2
		constterm =0;
		diffterm = 0;
		for x in range(0,len(cofficient)):
			constterm+= p2[x]*cofficient[x]
			diffterm+= diff[x]*cofficient[x]
		if((diffterm!=0).all()):
			t = ((-1)*constant - constterm)/diffterm
			pnt = t*(pt-p2)+p2
			pts1.append(pnt)
		else:
			pts1.append(pt)
	pts1 = np.array(pts1)
	#ax.scatter3D(pts1[:,0],pts1[:,1],pts1[:,2],color = "red")
	return pts1
	
def project_on_Sphere(point,cc,cr):
	"""
	Project points on SPhere
	"""
	projection = []
	for x in point:
		x = np.array(x)
		projection.append(((((x-cc)*cr)/math.dist(cc,x)))+cc)
	projection = np.array(projection)
	return projection
	
def generatesimplexfacets(poly,dimen):
	"""
	generate simplices of dimension dimen
	""" 
	tp = list(combinations(poly, dimen))
	listreturn = []
	for x in tp:
		listreturn.append(list(x))
	return listreturn

def hullfromtriangulation(simplices):
	"""
	Compute Hull From Triangulation
	"""
	hull= []
	for x in simplices:
		facets = generatesimplexfacets(x,len(x)-1)
		for f in facets:
			valid = True
			for y in simplices:
				if(valid==False):
					break
				if(x !=y).any():
					if(len(intersection(f,y))==len(f)):
						valid=False
						break
			if(valid==True):
				if(sorted(f) not in hull):
					hull.append(sorted(f))
	return hull
def computeChebshev_Sphere(points):
	"""
	Compute Chebshev Sphere of points
	"""
	polytophalfspaces = polyfun.compute_polytope_halfspaces(points)
	polyhalfspace = pc.Polytope(polytophalfspaces[0],polytophalfspaces[1])
	chebR = polyhalfspace.chebR
	chebXc = polyhalfspace.chebXc
	return chebXc,chebR

def plt_sphere(list_center, list_radius):
	"""
	Plot Spheres
	"""
	for c, r in zip(list_center, list_radius):
		#ax = fig.gca(projection='3d')
		# draw sphere
		u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
		x = r*np.cos(u)*np.sin(v)
		y = r*np.sin(u)*np.sin(v)
		z = r*np.cos(v)
		ax.plot_surface(x+c[0], y+c[1], z+c[2], color='red', alpha=0.2)

def checkplotsimplices(simplices,inputpoints):
	"""
	Plot Mesh Surface
	"""
	newlist = []
	i = 1
	for y in hull:
		pts = [inputpoints[x1] for x1 in y]
		ptrans = np.transpose(pts)
		verts = [list(zip(ptrans[0],ptrans[1],ptrans[2]))]
		ax.scatter3D(ptrans[0],ptrans[1],ptrans[2])
		ax.add_collection3d(Poly3DCollection(verts,facecolors='b', edgecolor = 'black', linewidths=1, alpha=.1))
def shift_Sphere_to_Orthogal_Axis(Points_on_Sphere,cc,cr):
	"""
	Orient Sphere to d-axis for convexization operations
	"""
	new_Sphere = []
	for x in Points_on_Sphere:
		newpoint = []
		for y in range(len(cc)-1):
			newpoint.append(x[y]-cc[y])
		newpoint.append(cr-x[y+1])
		new_Sphere.append(newpoint)
	return np.array(new_Sphere)

def closestAxis(point,diags):
	"""
	Find Closest Axis to Point
	"""
	i = 0
	assignedpartion = 0
	minvalue = 99999
	for B in diags:
		cosine =  np.arccos(np.dot(point,B)/(norm(point)*norm(B)))
		if(minvalue>cosine):
			minvalue = cosine
			assignedpartion = i
		i = i+1
	return assignedpartion
	
def rescalePoints(pts1,newpoints,diagnoal,chebR):
	"""
	Scale points of Steriographic Projection as if Peeling an Orange
	"""
	rescaledPoints = []
	radius = []
	for i in range(len(pts1)):
		axis = closestAxis(pts1[i],diagnoal)
		projectedPoint = [newpoints[i][j] for j in range(d-1)]
		cosine =  np.arccos(np.dot(projectedPoint,diagnoal[axis])/(norm(projectedPoint)*norm(diagnoal[axis])))
		prjectedRadius = norm(projectedPoint)
		arlength = (cosine*(2*np.pi*prjectedRadius))/360
		distance1 = norm(newpoints[i])
		center = [0 for j in range(d-1)]
		base = [0 for j in range(d)]
		center.append(chebR)
		center = np.array(center)
		base = np.array(base)
		theta = np.degrees(np.arccos(np.dot(base-center,newpoints[i]-center)/(norm(base-center)*norm(newpoints[i]-center))))
		distance = (theta*(2*np.pi*chebR))/360
		difference = distance - norm(newpoints[i])
		rescaledAngle = (prjectedRadius*cosine)/distance1
		pntsUI = []
		pntsUI.append([0 for j in range(d-1)])
		pntsUI.append(pts1[i])
		y = project_on_Sphere([diagnoal[axis]],[0 for j in range(d-1)],distance1)
		pntsUI.append(y[0])
		mean = (pntsUI[0]+pntsUI[1]+pntsUI[2])/3
		updatedpoints,pca = PCAprojectiondim(pntsUI,2)
		newupdatedpoints = []
		shift = updatedpoints[0]
		for x in updatedpoints:
			newupdatedpoints.append(x-updatedpoints[0])
		mat = rotation(np.array([1,0]),newupdatedpoints[2])
		mat2 = rotation(newupdatedpoints[2],np.array([1,0]))
		updatedpoints = []
		for x in newupdatedpoints:
			updatedpoints.append(mat.dot(x))
		negative = False
		if(updatedpoints[1][1]<0):
			negative = True
		x = distance*math.cos(rescaledAngle)
		y = distance*math.sin(rescaledAngle)
		if(negative):
			y = y*(-1)
		updatedpoints.append([x,y])
		newupdatedpoints = []
		for x in updatedpoints:
			newupdatedpoints.append(mat2.dot(x))
		updatedpoints = []
		for x in newupdatedpoints:
			updatedpoints.append(x+shift)
		data_original = np.dot(updatedpoints, pca.components_) +  mean # inverse_transform
		rescaledPoints.append(data_original[3])
		radius.append(difference)
	return np.array(rescaledPoints),np.array(radius)


def plotpolytop1(polyparts,d,inputpoints,name,radius):
	if(d==2):
		'''
		newpolyparts = set()
		edges = set()
		for x in polyparts:
			#print(x)
			#input()
			if(len(x)>5):
				for t in ConvexHull([inputpoints[i] for i in x],qhull_options="QJ").simplices:
					t = [x[i] for i in t]
					edges.add(tuple(sorted(t)))
					if(sorted(t) not in list(newpolyparts)):
						newpolyparts.add(tuple(sorted(t)))
		polyparts = newpolyparts
		'''
		hull = polyparts#ConvexHull(points).simplices
		fig, ax = plt.subplots()
	if(d==3):
		'''
		newpolyparts = []
		for x in polyparts:
			polysubparts = []
			projectedPoints = PCAprojection([inputpoints[i] for i in x])
			#for t in ConvexHull([inputpoints[i] for i in x],qhull_options="QJ").simplices:
			for t in Delaunay(projectedPoints).simplices:
				#print(t)
			#	print(Delaunay(projectedPoints).simplices)
				t = [x[i] for i in t]
			#	print(t)
			#	print(newpolyparts)
			#	input()
				if(sorted(t) not in polysubparts):
					polysubparts.append(sorted(t))
			newpolyparts.append(polysubparts)
		'''
		polyparts = polyparts
		hull = polyparts#ConvexHull(points).simplices
		fig = plt.figure()
		ax = Axes3D(fig, auto_add_to_figure=False)
		fig.add_axes(ax)
	cmap = plt.cm.get_cmap('hsv', len(hull))
	i = -1
	for y in hull:
		i = 	i+1
		if(d==3):
			if(len(y)>1):
				#fig = plt.figure()
				#ax = Axes3D(fig, auto_add_to_figure=False)
				#fig.add_axes(ax)
				np.random.seed()
				col = (np.random.rand(3,)+np.random.rand(3,)+np.random.rand(3,))/3
				for x in y:
					pts = [inputpoints[x1] for x1 in x]
					ptrans = np.transpose(pts)
					verts = [list(zip(ptrans[0],ptrans[1],ptrans[2]))]
					#ax.scatter3D(ptrans[0],ptrans[1],ptrans[2],s=5,color = 'red')
					ax.add_collection3d(Poly3DCollection(verts,facecolors=cmap(i),edgecolors = cmap(i), alpha=1))
				#plt.show()
		else:
			edges = [list(y)]
			for edge in edges:
				plt.plot([inputpoints[edge[0]][0],inputpoints[edge[1]][0]],[inputpoints[edge[0]][1],inputpoints[edge[1]][1]],color=cmap(i),linestyle='solid',linewidth=1)
			ptrans = np.transpose(inputpoints)
			patches = []
			for c,r in zip(inputpoints,radius):
				circle = Circle(c, r)
				patches.append(circle)
  
			# add these circles to a collection
			p = PatchCollection(patches)
			ax.add_collection(p)
			ax.scatter(ptrans[0],ptrans[1],s=5,color = 'red')	
	plt.savefig(name+".pdf", bbox_inches = 'tight',pad_inches = 0)
	plt.show()
	
def plotpolygon(polygons,dim,inputdata,name):
	fig, ax = plt.subplots(1,1)
	ax.set_xlim(-50,60)
	ax.set_ylim(-40,50)
	color1 = ['red','blue','green','yellow','cyan','pink']
	i=0
	for poly in polygons:
		if(len(poly)>10):
			simplices = Delaunay([inputdata[i] for i in poly]).simplices
			for s in simplices:
				polygon1 = Polygon([inputdata[poly[s[0]]],inputdata[poly[s[1]]],inputdata[poly[s[2]]]],color = color1[i%6],alpha=1)
				ax.add_patch(polygon1)
			i = i+1
	plt.savefig(name+".pdf", bbox_inches = 'tight',pad_inches = 0)
	plt.show()

def reduceDimension(points):
	"""
	Remove Last Coorindate as they are already projected on Plane
	"""
	pts = []
	for x in points:
		pts.append([x[i] for i in range(len(x)-1)])
	return np.array(pts)
	
def transformingConvexPolytopeForConvexDecomposition(dsphere):
	#******************Compute Chebshev Sphere *********************
	chebXc,chebR = computeChebshev_Sphere(dsphere)
	list_center = [(chebXc)]
	list_radius = [chebR]
	#*******************Project points on Sphere********************
	point_on_sphere = project_on_Sphere(dsphere,chebXc,chebR)

	#*******************Shift Sphere to orthogonal axis and touching the origin for Further Operations****************
	OrthoSphere = shift_Sphere_to_Orthogal_Axis(point_on_sphere,chebXc,chebR)
	#ax.scatter3D(OrtoSphere[:,0],OrtoSphere[:,1],OrtoSphere[:,2],s=10,color="blue")

	#*********** Choose two Projection Points across orthogonal axis with one being the origin and second being on the opposite side **************
	pp1 = np.array([0 for i in range(d)])
	pp2 = np.array([0 if i < (d-1) else chebR*2 for i in range(d)])
	pp1[0] = -0.001
	return OrthoSphere, pp1,pp2

def pointsToSimplicesEpsilon(OrtoSphere,center,d):
	
	newpoints = []
	diagnoal = simplexDiagonalWithCenter(0.0001,d)
	for pt in OrtoSphere:
		mat = rotation(pt-center,-center)
		for point in diagnoal:
			newpoints.append(mat.dot(point-center))
	newpoints = np.array(newpoints)
	return newpoints

def rotateSphere180deg(points,center):
	rotatedPoints = []
	pt = [i for i in center]
	pt[0] = 0.0001
	mat = rotation(-center,pt)
	for pt in points:
		rotatedPoints.append(mat.dot(pt)+(2*center))
	return np.array(rotatedPoints)

def equidistantFromObserver(points,r):
	azumutalFromOrigin = []
	for pt in points:
		length = 2*r*np.arctan(norm(pt)/(2*r))
		azumutalFromOrigin.append(((pt)*length)/norm(pt))
	return np.array(azumutalFromOrigin)
	
def LunicScaling(axis,azimPoint,origin,arcLen):
	pntsUI = np.array([origin,azimPoint,axis])
	mean = (pntsUI[0]+pntsUI[1]+pntsUI[2])/3
	updatedpoints, pca = PCAprojectiondim(pntsUI,2)
	newupdatedpoints = []
	shift = updatedpoints[0]
	for x in updatedpoints:
		newupdatedpoints.append(x-updatedpoints[0])
	mat = rotation(np.array([1,0]),newupdatedpoints[2])
	mat2 = rotation(newupdatedpoints[2],np.array([1,0]))
	updatedpoints = []
	for x in newupdatedpoints:
		updatedpoints.append(mat.dot(x))
	negative = False
	if(updatedpoints[1][1]<0):
		negative = True
	distance = norm(np.array(azimPoint))
	rescaledAngle = arcLen/distance
	x = distance*math.cos(rescaledAngle)
	y = distance*math.sin(rescaledAngle)
	if(negative):
		y = y*(-1)
	updatedpoints.append([x,y])
	newupdatedpoints = []
	for x in updatedpoints:
		newupdatedpoints.append(mat2.dot(x))
	updatedpoints = []
	for x in newupdatedpoints:
		updatedpoints.append(x+shift)
	data_original = np.dot(updatedpoints, pca.components_) +  mean # inverse_transform
	return data_original[3]
		
	
def equidistantFromObserverAxis(points,originalPoints,axis):
	lunePoints = []	
	origin = np.array([0 for i in range(len(points[0]))])
	for i in range(len(points)):
		projectedPoint = np.array([originalPoints[i][j] for j in range(len(points[0]))])
		radius = norm(projectedPoint)
		theta = np.arccos(np.dot(axis,projectedPoint)/(norm(axis)*norm(projectedPoint)))
		arcLen = radius*theta
		projectionAxis = (axis*norm(points[i]))/norm(axis)
		recaledPoint = LunicScaling(projectionAxis,points[i],origin,arcLen)
		lunePoints.append(recaledPoint)
	return np.array(lunePoints)
		
def cylindricalProjection(points,center,t):
	#Project on Cylinder
	radius = norm(center)
	origin = np.array([0 for i in center])
	orthogonalPoint = np.array([i for i in center])
	orthogonalPoint[0] = radius
	cylinricalprojection= []
	vectorOnCilindricalSurface = []
	cylinricalprojectiononPlane = []
	proj = []
	vectorOnplane = []
	for i in range(len(points)):
		projectedPoint = [points[i][j] for j in range(len(points[0])-1)]
		projectedPoint.append(radius)
		projectedPoint = np.array(projectedPoint)
		proj.append(projectedPoint)
		projectionOnSphere = (((projectedPoint-center)*radius)/math.dist(center,projectedPoint))+center
		theta = np.arccos(np.dot(points[i]-center,projectedPoint-center)/(norm(points[i]-center)*norm(projectedPoint-center)))
		radialAngle =  np.arccos(np.dot(orthogonalPoint-center,projectedPoint-center)/(norm(orthogonalPoint-center)*norm(projectedPoint-center)))
		radialDist = radius*radialAngle
		if(t==1):
			radialDist = np.pi*radius - radialDist
		r = radius/np.cos(theta)
		pointOnCylinder  = (((points[i]-center)*r)/math.dist(center,points[i]))+center
		vectorOnCilindricalSurface.append(pointOnCylinder-projectionOnSphere)
		cylinricalprojection.append(pointOnCylinder)
		if(points[i][0]<0):
			radialDist = -radialDist
		cylinricalprojectiononPlane.append(radialDist)
	projectedPoints,pca = PCAprojectiondim(vectorOnCilindricalSurface,len(vectorOnCilindricalSurface[0])-2)
	for x,r in zip(projectedPoints,cylinricalprojectiononPlane):
		pt = []
		pt.append(r)
		for coord in x:
			pt.append(coord)
		vectorOnplane.append(pt)
	return np.array(proj),np.array(cylinricalprojection),np.array(vectorOnplane)
	
def pruneMaximalParts(convexFaces,convexpartssorted,delaunayparts,weights,dim):
	""" Extract Maximal Components """
	for i in range(len(delaunayparts)):
		polytope = delaunayparts[i]
		present = False
		done = False
		for j in range(len(convexFaces)):
			candidatePolytope = convexFaces[j]
			if(len(polytope)>0):
				for x in polytope:
					stop = False
					for y in candidatePolytope:
						if x==y:
							present = True
							if len(polytope)>len(candidatePolytope):
								convexFaces[j] = polytope
								#done = True
							#stop = True
					#if(stop):
						#break
			else:
				present = True
				break
			if(done):
				break
			if(present):
				break	
		if(not present):
			convexFaces.append(polytope)
	return convexFaces
			
def generateConvexFaces(convexFaces,hull,newpoints,pp1,pp2):
	d = len(newpoints[0])-1
	radius = [0 for i in range(len(newpoints))]
	'''
	fig = plt.figure()
	ax = Axes3D(fig, auto_add_to_figure=False)
	fig.add_axes(ax)	
	ax.scatter3D(newpoints[:,0],newpoints[:,1],newpoints[:,2],color = "blue")
	plt.show()
	'''
	#*****************First Step Compute StereoGraphic Projection********************
	sterographicProjection = reduceDimension(projectonplane(pp1,pp2,newpoints))
	convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,d,sterographicProjection)
	convexFaces = pruneMaximalParts(convexFaces,convexpartssorted,delaunayparts,weights,d)
	#plotpolytop1(delaunayparts,len(sterographicProjection[0])+1,newpoints,"name1",radius)
	#*****************Second Step Compute Azimutal Projection*********************
	azimuthalProjection = equidistantFromObserver(sterographicProjection,norm(pp2/2))
	convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,d,azimuthalProjection)
	convexFaces = pruneMaximalParts(convexFaces,convexpartssorted,delaunayparts,weights,d)
	#plotpolytop1(delaunayparts,len(sterographicProjection[0])+1,newpoints,"name1",radius)
	#*****************Projection on Lune with central Axis ***************************
	diagonal = simplexDiagonals(1,d+1)
	i=0
	for axis in diagonal:
		luneProjection = equidistantFromObserverAxis(azimuthalProjection,newpoints,axis)
		convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,d,luneProjection)
		convexFaces = pruneMaximalParts(convexFaces,convexpartssorted,delaunayparts,weights,d)
		#plotpolytop1(delaunayparts,len(sterographicProjection[0])+1,newpoints,"name1",radius)
		i = i+1
	#*************** Finally Compute the Cylindrical projection on to Plane***************
	proj,c3d,cProjection = cylindricalProjection(newpoints,pp2/2,0)
	convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,d,cProjection)
	convexFaces = pruneMaximalParts(convexFaces,convexpartssorted,delaunayparts,weights,d)
	proj,c3d,cProjection = cylindricalProjection(newpoints,pp2/2,1)
	convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,d,cProjection)
	convexFaces = pruneMaximalParts(convexFaces,convexpartssorted,delaunayparts,weights,d)
	
	#plotpolytop1(delaunayparts,len(sterographicProjection[0])+1,newpoints,"name1",radius)

	'''
	plotpolytop1(convexpartssorted,len(cProjection[0]),cProjection,"name4"+str(5),radius)
	fig = plt.figure()
	ax = Axes3D(fig, auto_add_to_figure=False)
	fig.add_axes(ax)	
	ax.scatter3D(newpoints[:,0],newpoints[:,1],newpoints[:,2],color = "blue")
	ax.scatter3D(c3d[:,0],c3d[:,1],c3d[:,2],color = "green")
	ax.scatter3D(proj[:,0],proj[:,1],proj[:,2],color = "cyan")
	plt.show()
	fig, ax = plt.subplots()
	ax.scatter(cProjection[:,0],cProjection[:,1],color = "red")
	plt.show()
	'''
	return convexFaces
def getUniqueVertices(face):
	uV = set()
	for x in face:
		for y in x:
			uV.add(y)
	return uV
def informedConvexization(convexFaces,hull,newpoints,pp1,pp2):
	i =0
	d = len(newpoints[0])-1
	radius = [0 for i in range(len(newpoints))]
	updatedConvexFaces = []
	#plotpolytop1(convexFaces,len(newpoints[0]),newpoints,"name1",radius)
	for face in convexFaces:
		if len(face)>2:
			uniqueVertices = getUniqueVertices(face)
			centroid = np.array([0 for i in range(len(newpoints[0]))])
			for x in uniqueVertices:
				centroid = centroid + newpoints[x]
			centroid = centroid/len(uniqueVertices)
			pp1 = (((centroid-pp2/2)*norm(pp2/2))/math.dist(pp2/2,centroid))+(pp2/2)
			pp2 = pp2 + pp1
			sterographicProjection = reduceDimension(projectonplane(pp1,pp2,newpoints))
			azimuthalProjection = equidistantFromObserver(sterographicProjection,norm(pp2/2))
			convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,d,azimuthalProjection)
			updatedConvexFaces = pruneMaximalParts(updatedConvexFaces,convexpartssorted,delaunayparts,weights,d)
			#plotpolytop1(updatedConvexFaces,len(sterographicProjection[0])+1,newpoints,"name1",radius)
		i = i +1
	return updatedConvexFaces

#**************** Input******************

#dsphere = np.array(fiblat.sphere_lattice(3,100))
#dsphere = tadasets.dsphere(n=100, d=2, r=1)
d=3
#dsphere = np.loadtxt("input3dim4size100noise0",delimiter = ",")
#dsphere = np.loadtxt("Tetrahedron.txt",delimiter = ",")
dsphere = np.loadtxt("input3dim3size150noise1",delimiter = ",")
computedistancematrix(dsphere)

#******************Compute Delaunay Triangulation***************
triangulation =  Delaunay(dsphere).simplices

#******************Compute Hull from Delaunay Triangulation***********
hull = hullfromtriangulation(triangulation)

newpoints, pp1, pp2 = transformingConvexPolytopeForConvexDecomposition(dsphere)

convexFaces = []
convexFaces = generateConvexFaces(convexFaces,hull,newpoints,pp1,pp2)
rotatedPoints = rotateSphere180deg(newpoints,pp2/2)
convexFaces = generateConvexFaces(convexFaces,hull,rotatedPoints,pp1,pp2)

convexFaces = informedConvexization(convexFaces,hull,newpoints,pp1, pp2)

count = 0
for x in convexFaces:
	if len(x) > 1:
		count = count +1
print(count)
print(len(convexFaces))

col = ["red","blue","green","cyan"]	
i = 0
radius = [0 for i in range(len(dsphere	))]
plotpolytop1(convexFaces,len(dsphere[0]),dsphere,"plotsphere"+str(5),radius)

fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)

for y in convexFaces:
	for x in y:
	#ax.scatter3D(dsphere[:,0],dsphere[:,1],dsphere[:,2],color = col[3])
		points = np.array([dsphere[i] for i in x])
		ax.scatter3D(points[:,0],points[:,1],points[:,2],color = col[i%4],s=15)
	i = i+1
plt.show()
fig, ax = plt.subplots()
cProjection = PCAprojection(points)
ax.scatter(cProjection[:,0],cProjection[:,1],color = "red")
plt.show()

#plt.show()
#pts2 = reduceDimension(projectonplane(pp1,pp2,newpoints))
pts1 = reduceDimension(projectonplane(pp1,pp2,newpoints))
convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,len(pts1[0]),pts1)
radius = [0 for i in range(len(pts1))]
plotpolytop1(convexpartssorted,len(pts1[0]),pts1,"name3",radius)
print(convexpartssorted[0])
print(convexpartssorted[1])
print(convexpartssorted)
print(delaunayparts)
print(weights)
#print(pp2)
#print(pp1)
diagnoal = simplexDiagonals(1,d)
print(diagnoal)
newvector = np.array([0 for i in range(d-1)])
for x in range(d-1):
	newvector = newvector + diagnoal[x]
newvector = newvector/(d-1)
print(newvector)
print(diagnoal[0])
#input()
mat = rotation(newvector,diagnoal[0])
diagonal2 = []
for x in diagnoal:
	diagonal2.append(mat.dot(x))

scaledPoints,radius1 = rescalePoints(pts1,newpoints,diagnoal,chebR)
scaledPoints2,radius2 = rescalePoints(pts1,newpoints,diagonal2,chebR)
convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,len(scaledPoints[0]),scaledPoints)
plotpolytop1(convexpartssorted,len(pts1[0]),scaledPoints,"name1",radius1)
convexpartssorted,delaunayparts,weights = iterativeconvexization(hull,len(scaledPoints2[0]),scaledPoints2)
plotpolytop1(convexpartssorted,len(pts1[0]),scaledPoints2,"name2",radius2)

plt.show()
phi = np.linspace(0, 2*np.pi, len(dsphere))
rgb_cycle = (np.stack((np.cos(phi          ), # Three sinusoids,
                       np.cos(phi+2*np.pi/3), # 120° phase shifted,
                       np.cos(phi-2*np.pi/3)
                      )).T # Shape = (60,3)
             + 1)*0.5                         # scaled to [0,1].
             
fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)	
pts = pts1
#ax.scatter3D(pts[:,0],pts[:,1],pts[:,2],s=10,color="green")
ax.scatter(pts[:,0],pts[:,1],color='green',s=10)
plt.show()
fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)	
#ax.scatter3D(scaledPoints[:,0],scaledPoints[:,1],scaledPoints[:,2],s=10,color=rgb_cycle)
ax.scatter(scaledPoints[:,0],scaledPoints[:,1],color=rgb_cycle,s=5)
plt.show()
fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)	
#ax.scatter3D(scaledPoints2[:,0],scaledPoints2[:,1],scaledPoints2[:,2],s=10,color=rgb_cycle)
ax.scatter(scaledPoints2[:,0],scaledPoints2[:,1],color=rgb_cycle,s=5)
plt.show()
# Change major ticks to show every 0.2.
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(5))

# Change minor ticks to show every 0.05. (0.2/4 = 0.05)
ax.xaxis.set_minor_locator(AutoMinorLocator(0.5))
ax.yaxis.set_minor_locator(AutoMinorLocator(0.5))

# Turn grid on for both major and minor ticks and style minor slightly
# differently.
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
figure = plt.gcf() # get current figure
figure.set_size_inches(20, 20)
plt.xticks(rotation=90, ha='right')
#patches = []
#for i in range(0,len(pts1),4):
#	center, radius  = define_circle(pts1[i],pts1[i+1],pts1[i+2])
#	if center is not None:
#		circle = Circle(center, radius)
#		patches.append(circle)
  
# add these circles to a collection
#p = PatchCollection(patches, color = "green", alpha = 0.1)
#ax.add_collection(p)
#ax.set_xlim(-50,650)
#ax.set_ylim(-450,50)
plt.savefig("outputPIpolytopalproj2.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()
