
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
import tadasets
from mpl_toolkits.mplot3d import Axes3D

dimension = 2
points = tadasets.dsphere(n=1000, d=dimension-1, r=1, noise=0.15)

#print(np.array(points))
#points1 = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1], [0.5,1.5],[0.5,-0.5],[2,2],[2.5,1],[3,-1],[1.5,-0.5],[1.7,-0.3],[1.9,0.2],[2.3,0.3],[2.4,-1.5],[2.3,1.2],[0.5,0.5],[2.2,5.6],[3.7,1.3]])
#plt.scatter(points[:,0],points[:,1])
#plt.show()
tri = Delaunay(points)
print(len(tri.simplices))
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))
 
def union(lst1,lst2):
	return list(set().union(lst1 , lst2))

def neighbours(polytop,simplices):
	polytop = list(polytop)
	adjacency = []
	for y in simplices:
		y = list(y)
		if(polytop != y):
			if(len(intersection(polytop,y))==dimension):
				adjacency.append(y)
	return adjacency
#xt = int(math.sqrt(len(simplices)))+1
#fig, axs = plt.subplots(xt,xt)

def mergeneighbors(polytop,simplices):
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


#for in adjacency
def optimalconvexization(simplices):
	convexparts = []
	maxdist = []
	averagedist = []
	sizecd = []
	for i in range(0,len(simplices)):
		#axs[int(i/xt)][i%xt].scatter(points[:,0],points[:,1])
		x= list(simplices[i])
		distance =0
		maxval = 0
		while True:
			x1 = list(mergeneighbors(x,simplices))
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


valid = []
def iterativeconvexization(simplices):
	convexpartsunion = []
	simplices = list(simplices)
	remainingsimplices = simplices
	premaining = 0
	while(True):
		df = optimalconvexization(remainingsimplices)
		i=0
		pointsaddressed = []
		for x in df['convexpart'].tolist():
			check = 1
			for t in convexpartsunion:
				if(len(intersection(t,x)) > dimension):
					check = 0
			if(check==1):
				valid.append(i)
				hull = ConvexHull([points[i] for i in x])
				#for p in hull.simplices:
				#	axs[int(i/xt)][i%xt].plot([points[x[p[0]]][0],points[x[p[1]]][0]],[points[x[p[0]]][1],points[x[p[1]]][1]])
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
		print(remaining)
		if(remaining==premaining):
			if(remaining !=0):
				for r in remainingsimplices:
					convexpartsunion.append(r)
			break
		premaining = remaining
	return convexpartsunion

convexpartsunion = iterativeconvexization(tri.simplices)

'''
i=0
remaining_edges = []
for x in df['convexpart'].tolist():
	if i in valid:
		pointscov = []
		pairvoc = []
		for y in convexpartsunion:
			if(x != y):
				coverededge = intersection(x,y)	
				pointscov = union(pointscov,coverededge)
				pairvoc.append(coverededge)
		hull = ConvexHull([points[i] for i in x])
		pt = (hull.simplices).tolist()
		pt = [sorted([x[t] for t in y]) for y in pt]
		pairvoc = [sorted(y) for y in pairvoc]
		pairvoc = [i for i in pt if i not in pairvoc]
		hull = ConvexHull(points)
		pt = (hull.simplices).tolist()
		pt = [sorted(y) for y in pt]
		pairvoc = [i for i in pairvoc if i not in pt]
		remaining_edges.append([x for x in pairvoc])
	i =i+1 	
j=15
print(remaining_edges)
newlist = []
reminingconvexparts = []
i=0
for x in remaining_edges:
	for y in remaining_edges[i+1:]:
		if(x!=y):
			for t in x:
				for p in y:
					if(len(intersection(t,p))>0):
						 newlist.append(t)
						 newlist.append(p)
						 reminingconvexparts.append(union(t,p))
	i =i+1
print(newlist)
ts =[]
print(reminingconvexparts)
for x in remaining_edges:
	for y in x:
		ts.append(y)

L = [i for i in ts if i not in newlist]

print(L)
remlist = []
for x in reminingconvexparts:
	hull = ConvexHull([points[i] for i in x])
	pt = (hull.simplices).tolist()
	pt = [sorted([x[t] for t in y]) for y in pt]
	print(pt)
	redge = [i for i in pt if i not in newlist]
	remlist.append(redge)

print(remlist)

	
print(newlist)
for x in L:
	for y in remlist:
		for tr in y:
			if(len(intersection(x,tr))>0):
				reminingconvexparts.append(union(x,tr))

#print(reminingconvexparts)
#print(convexpartsunion)
for t in reminingconvexparts:
	for y in convexpartsunion:
		if(len(intersection(t,y))<3):
			convexpartsunion.append(t)
'''
	

#X=np.array(points)
#fig = plt.figure(figsize=(10,7))
#ax=Axes3D(fig)
#ax.plot([i,i],[j,j],[k,h],color = 'g')
#plt.show()
#X1 = np.transpose(X)[0]
#Y1 = np.transpose(X)[1]
#Z1 = np.transpose(X)[2]


for x in convexpartsunion:
	print(len(x))
	hull = ConvexHull([points[i] for i in x])
	for p in hull.simplices:
		for t in range(0,len(p)):
			for e in range(t+1,len(p)): 
				plt.plot([points[x[p[0]]][0],points[x[p[1]]][0]],[points[x[p[0]]][1],points[x[p[1]]][1]])
				#ax.plot([X1[x[p[t]]],X1[x[p[e]]]], [Y1[x[p[t]]],Y1[x[p[e]]]], [Z1[x[p[t]]],Z1[x[p[e]]]])
				 
plt.scatter(points[:,0],points[:,1],s=10,color='black')
#ax.scatter(np.transpose(X)[0], np.transpose(X)[1], np.transpose(X)[2])

plt.show()


'''
for x in adjacency:
	print(x)
	
for x in adjacency:
	adj =[]
	for y in x:
		hull = ConvexHull([points[i] for i in y])
		uniquedict = {}
		for sublist in hull.simplices:
			for item in sublist:
				if item in uniquedict.keys():
					uniquedict[item] += 1
				else:
					uniquedict[item] = 1
		if(len(list(uniquedict.keys())) == len(y)):
			print("True")
			mhull = []
			for i in hull.simplices:
				tmp = []
				for j in i:
					tmp.append(y[j])
				mhull.append(tmp)
			adj.append(mhull)
		else:
			print("False")
			adjacency[index].remove(y)
	index= index+1
	adjacency2.append(adj)
		  
for x in adjacency2:
	print(x)

for x in adjacency2:
	for y in x:
		for p in y:
			plt.plot([points[p[0]][0],points[p[1]][0]],[points[p[0]][1],points[p[1]][1]])
plt.show()
index=0
adjacency4 = []
for x in adjacency:
	adj = []
	tp =0
	for y in x:
		tp=tp+1
		for z in x[tp:]:
			ut =union(y,z)
			hull = ConvexHull([points[i] for i in ut])
			uniquedict = {}
			for sublist in hull.simplices:
				for item in sublist:
					if item in uniquedict.keys():
						uniquedict[item] += 1
					else:
						uniquedict[item] = 1
			if(len(list(uniquedict.keys())) == len(ut)):
				adjacency4.append(list(uniquedict.keys()))
				print("True")
				mhull = []
				for i in hull.simplices:
					tmp = []
					for j in i:
						tmp.append(ut[j])
					mhull.append(tmp)
				adj.append(mhull)
			else:
				print("False")  
	index= index+1
	adjacency3.append(adj)			

print(len(adjacency3))		
for x in [adjacency3[i] for i in [0,1,2]]:
	for y in x:
		for p in y:
			plt.scatter(points[:,0],points[:,1])
			plt.plot([points[p[0]][0],points[p[1]][0]],[points[p[0]][1],points[p[1]][1]])
	#input("x")		
	#plt.show()
plt.show()
	
adjacency5 =[]
index=0
for x in [1]:
	adj = []
	tp =0
	for y in adjacency4:
		tp=tp+1
		for z in adjacency4:
			ut =union(y,z)
			if(len(intersection(y,z))==2):
				hull = ConvexHull([points[i] for i in ut])
				uniquedict = {}
				for sublist in hull.simplices:
					for item in sublist:
						if item in uniquedict.keys():
							uniquedict[item] += 1
						else:
							uniquedict[item] = 1
				if(len(list(uniquedict.keys())) == len(ut)):
					print("True")
					mhull = []
					for i in hull.simplices:
						tmp = []
						for j in i:
							tmp.append(ut[j])
						mhull.append(tmp)
					adj.append(mhull)
				else:
					print("False")  
	index= index+1
	adjacency5.append(adj)			

print(adjacency5)
for x in [1]:
	for y in adjacency5:
		for p in y:
			plt.scatter(points[:,0],points[:,1])
			plt.plot([points[p[0]][0],points[p[1]][0]],[points[p[0]][1],points[p[1]][1]])
		input("x")		
		plt.show()
'''
