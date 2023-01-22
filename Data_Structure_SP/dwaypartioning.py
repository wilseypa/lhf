import itertools
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D      
import matplotlib
import tadasets
import fiblat

print("Enter Number of Points")
pts = int(input())
print("Enter Dimension of Points")
dim = int(input())

from pca import pca

#k = np.random.random ([pts,dim])
#np.savetxt("inputfile.txt",k)
k = np.loadtxt("inputfile.txt")
#print(k)
#input()
#k = tadasets.dsphere(n=pts,d=dim-1, r=5,noise = 0)
#k =  np.array(fiblat.sphere_lattice(dim,pts))

points = k

maxaxislen = 10

'''
from scipy.spatial import Delaunay
tri = Delaunay(points)
plt.triplot(points[:,0], points[:,1], tri.simplices,lw=0.5)
plt.plot(points[:,0], points[:,1], 'o',markersize=1)
plt.savefig("Delaunay.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()
'''

hyplane = [[0 if k!= i else 20 for k in range(dim+1)] for i in range(dim+1)]
model = pca(n_components=dim)
coords = model.fit_transform(hyplane)
coords = np.array([coords["PC"].iloc[:, i] for i in range(dim)]).T

class dwayTreeNode:
	def __init__(self, centroid,data):
		self.parent = None
		self.centroid = centroid
		self.data = data
		self.children = []


def inorder(node):
        if node == None:
            return

        for child in node.children:
            inorder(child)
        datapoints.append(node.centroid)

def buildDwayTree(data):
	data = np.array(data)
	if(data.size<=0):
		return None
	centroid = data.mean(axis=0)
	nNode = dwayTreeNode(centroid,[])
	if (nNode == None):
		print("Memory error")
		return None
	
	partitions = []
	for A in data:
		A_vec =A-centroid
		i = 0
		assignedpartion = 0
		maxvalue = 0
		for B in coords:
			cosine =  np.arccos(np.dot(A_vec,B)/(norm(A_vec)*norm(B)))*90
			if(maxvalue<cosine):
				maxvalue = cosine
				assignedpartion = i
			i = i+1
		partitions.append(assignedpartion)
	
	dwaypartition = [[] for i in range(dim+1)]
	i = 0
	for x in partitions:
		dwaypartition[x].append(list(data[i]))
		i = i+1
	for part in dwaypartition:
		if(len(part)>dim+1):
			child = buildDwayTree(part)
			child.parent = nNode
			nNode.children.append(child)
		if(len(part)<=dim+1):
			data1 = np.array(part)
			if(data1.size>0):
				centroid1 = data1.mean(axis=0)
			else:
				centroid1 = None
			child = dwayTreeNode(centroid1,part)
			child.parent = nNode
			nNode.children.append(child)
	return nNode

def get_neighbor_of_greater_or_equal_size(node, direction):   
	if node.parent is None:
		return None
	if node.parent.children[direction] != node: 
		return node.parent.children[direction]
		
	otherneighbornode = get_neighbor_of_greater_or_equal_size(node.parent,direction)
def	is_leaf(node):
	if node is None:
		return True
	if(len(node.data)<=dim+1):
		return True
	else:
		return False
	
def find_neighbors_of_smaller_size(node, neighbor, direction):   
	candidates = [] if neighbor is None else neighbor
	neighbors = []
	while len(candidates) > 0:
		if is_leaf(candidates[0]):
			neighbors.append(candidates[0])
		else:
			for i in range(1,dim+1):
				if(candidates[0].children[direction+i%(dim+1)] is not None):
					candidates.append(candidates[0].children[direction+1%dim])
		candidates.remove(candidates[0])
	return neighbors
	
def get_neighbors(node, direction):   
	neighbors = []
	firstneigh = get_neighbor_of_greater_or_equal_size(node,direction)
	if firstneigh is not None:
		neighbors.append(firstneigh) 		#node not a parents i child, parents i child is always a i direction neighbour
	secondneigh = get_neighbor_of_greater_or_equal_size(node.parent,direction)	#parent->parent until node is not a i child --> parent i child is a neighbor	
	if secondneigh is not None:
		neighbors.append(secondneigh) 		#node not a parents i child, parents i child is always a i direction neighbour
	finalneighbors = find_neighbors_of_smaller_size(node,neighbors, direction)
	return finalneighbors
    	
root = buildDwayTree(points)
node = root.children[2].children[1].children[1]
neighbor = []
print(node.centroid,"sdf")
for i in range(dim+1):
	for y in get_neighbors(node,i):
		neighbor.append(y)
for neigh in neighbor:
	print(neigh.centroid)

datapoints = []
inorder(root)


for x in datapoints:
	print(x)
	print("******")
plt.plot(points[:,0], points[:,1], 'o',markersize=1)
plt.show()
'''
#fig, ax = plt.subplots()
#fig.canvas.draw()  # if running all the code in the same cell, this is required for it to work, not sure why
xcoords = []
ycoords = []
zcoords = []
colorlist = []
count = 1
si = []
lines = []
for x in datapoints:
	if(type(x[0])==type([])):
		for y in x:
			y = list(y)
			xcoords.append(y[0])
			ycoords.append(y[1])
			zcoords.append(y[2])
			colorlist.append(count)
			si.append(10)
	else:
		x = list(x)
		xcoords.append(x[0])
		ycoords.append(x[1])
		zcoords.append(x[2])
		colorlist.append(0)
		#lines.append(Line2D([x[0],x[0]+coords[0][0]], [x[1],x[1]+coords[0][1]],linewidth= 0.1))
		#lines.append(Line2D([x[0],x[0]+coords[1][0]], [x[1],x[1]+coords[1][1]],linewidth= 0.1))
		#lines.append(Line2D([x[0],x[0]+coords[2][0]], [x[1],x[1]+coords[2][1]],linewidth= 0.1))
		si.append(0)
		count = count+1
		
colorlist = list(np.array(colorlist)/len(colorlist))
#ax.scatter(xcoords,ycoords,cmap='inferno', c=colorlist,s=si)
#for l in lines:
#	ax.add_line(l)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xcoords,ycoords,zcoords,cmap='inferno', c=colorlist,s=si)
plt.savefig("dway.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()
'''
'''
c=["red","blue","green"]
i = 0
fig = plt.figure()
ax = fig.add_subplot(111)
for x in dwaypartition:
	if(len(x)>0):
		x = np.array(x)
		ax.scatter(x[:,0], x[:,1], s=1,color = c[i])
	i = i+1
line1 = Line2D([0+centroid[0],coords[0][0]+centroid[1]], [0+centroid[0],coords[0][1]+centroid[1]])
line2 = Line2D([0+centroid[0],coords[1][0]+centroid[1]], [0+centroid[0],coords[1][1]+centroid[1]])
line3 = Line2D([0+centroid[0],coords[2][0]+centroid[1]], [0+centroid[0],coords[2][1]+centroid[1]])
ax.add_line(line1)
ax.add_line(line2)
ax.add_line(line3)
'''
	
	
