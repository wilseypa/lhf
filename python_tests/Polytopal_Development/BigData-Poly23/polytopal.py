import tadasets
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
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
	for x in point_on_sphere:
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
	
def rescalePoints(pts1,newpoints,diagnoal):
	"""
	Scale points of Steriographic Projection as if Peeling an Orange
	"""
	rescaledPoints = []
	for i in range(len(pts1)):
		axis = closestAxis(pts1[i],diagnoal)
		projectedPoint = [newpoints[i][j] for j in range(d-1)]
		cosine =  np.arccos(np.dot(projectedPoint,diagnoal[axis])/(norm(projectedPoint)*norm(diagnoal[axis])))
		prjectedRadius = norm(projectedPoint)
		distance = norm(pts1[i])
		rescaledAngle = (prjectedRadius*cosine)/distance
		pntsUI = []
		pntsUI.append([0 for j in range(d-1)])
		pntsUI.append(pts1[i])
		y = project_on_Sphere([diagnoal[axis]],[0 for j in range(d-1)],distance)
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
	return np.array(rescaledPoints)
	
def reduceDimension(points):
	"""
	Remove Last Coorindate as they are already projected on Plane
	"""
	pts = []
	for x in points:
		pts.append([x[i] for i in range(len(x)-1)])
	return np.array(pts)
	
fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)	

dsphere = np.array(fiblat.sphere_lattice(3,200))
#dsphere = tadasets.dsphere(n=200, d=2, r=1)
d=3
#dsphere = np.loadtxt("Tetrahedron.txt",delimiter = ",")
dsphere = np.loadtxt("input0dim3size500noise1.txt",delimiter = ",")
triangulation =  Delaunay(dsphere).simplices
hull = hullfromtriangulation(triangulation)
chebXc,chebR = computeChebshev_Sphere(dsphere)
list_center = [(chebXc)]
list_radius = [chebR]
#checkplotsimplices(hull,dsphere)
point_on_sphere = project_on_Sphere(dsphere,chebXc,chebR)
#ax.scatter3D(point_on_sphere[:,0],point_on_sphere[:,1],point_on_sphere[:,2],s=10,color="red")
OrtoSphere = shift_Sphere_to_Orthogal_Axis(point_on_sphere,chebXc,chebR)
ax.scatter3D(OrtoSphere[:,0],OrtoSphere[:,1],OrtoSphere[:,2],s=10,color="blue")

#ProjectionPoints
pp1 = np.array([0 for i in range(d)])
pp2 = np.array([0 if i < (d-1) else chebR*2 for i in range(d)])
pp1[0] = -0.001
'''
newpoints = []
diagnoal = simplexDiagonalWithCenter(0.0001,d)
for pt in OrtoSphere:
	mat = rotation(pt-pp2/2,-pp2/2)
	for point in diagnoal:
		newpoints.append((mat.dot(point-pp2/2))+pp2/2)
newpoints = np.array(newpoints)
'''
newpoints = OrtoSphere
ax.scatter3D(newpoints[:,0],newpoints[:,1],newpoints[:,2],s=10,color="green")
#pts2 = reduceDimension(projectonplane(pp1,pp2,newpoints))
pts1 = reduceDimension(projectonplane(pp1,pp2,newpoints))

plt.show()	
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

scaledPoints = rescalePoints(pts1,newpoints,diagnoal)
scaledPoints2 = rescalePoints(pts1,newpoints,diagonal2)

plt.show()
phi = np.linspace(0, 2*np.pi, 456)
rgb_cycle = (np.stack((np.cos(phi          ), # Three sinusoids,
                       np.cos(phi+2*np.pi/3), # 120° phase shifted,
                       np.cos(phi-2*np.pi/3)
                      )).T # Shape = (60,3)
             + 1)*0.5                         # scaled to [0,1].
             
fig, ax = plt.subplots()
pts = pts1
ax.scatter(pts[:,0],pts[:,1],color='green',s=10)
ax.scatter(scaledPoints[:,0],scaledPoints[:,1],color=rgb_cycle,s=5)
ax.scatter(scaledPoints2[:,0],scaledPoints2[:,1],color=rgb_cycle,s=5)

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
