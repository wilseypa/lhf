import  tadasets
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn.decomposition import PCA
import random


def PCAprojection(pts):
	pca = PCA(n_components=len(pts[0])-1).fit(pts)
	pca_x = pca.transform(pts)
	return pca_x
	

datasize  = 10
dim = 3
unitsphere = np.array([np.array([ 0.70013219, -0.6740024 , -0.23566008]), np.array([ 0.89640578, -0.4363997 ,  0.07753701]), np.array([ 0.41702517, -0.85432765,  0.31018425]), np.array([0.19171654,0.31713976, 0.92879876]), np.array([-0.33456917,  0.82684448, -0.45209698]), np.array([ 0.35478838, -0.9346821 , -0.02223899]), np.array([-0.36459425,  0.82332563, -0.43497809]), np.array([-0.72403641, -0.08169234,  0.68490703]), np.array([-0.64549831,  0.5052651 , -0.572747  ]), np.array([-0.47213208,  0.67731401,  0.56421364])])  #tadasets.dsphere(n=datasize,d=dim-1,r=1,noise=0)
centroid = sum(unitsphere)/len(unitsphere)
#print(centroid)
newunitsphere = [x-centroid for x in unitsphere]
#ax = plt.axes(projection='3d')
'''
zdata = [item[0] for item in unitsphere]
xdata = [item[1] for item in unitsphere]
ydata = [item[2] for item in unitsphere]
ax.scatter3D(xdata, ydata, zdata,color = 'blue')
ax.scatter3D([centroid[0]],[centroid[1]],[centroid[2]],color = 'blue')

zdata = [item[0] for item in newunitsphere]
xdata = [item[1] for item in newunitsphere]
ydata = [item[2] for item in newunitsphere]
ax.scatter3D(xdata, ydata, zdata,color = 'red')
centroid = [0,0,0]
ax.scatter3D([centroid[0]],[centroid[1]],[centroid[2]],color = 'red')
'''
def simplexpoints(c):
	newpoints = []
	cnt = 0
	for x in range(0,c):
		coords = []
		for y in range(0,c):
			if(y==cnt):
				coords.append(1)
			else:
				coords.append(0)
		newpoints.append(coords)
		cnt = cnt+1
	dimlesspoints = PCAprojection(newpoints)
	updatedsphere = []
	for x in dimlesspoints:
		updatedsphere.append((0.01*x)/(math.dist([0,0],x)))
	return(updatedsphere)
def getsimplexsphere(points):
	c = len(points[0])
	pca = PCA(n_components=c-1).fit(points)
	pca_x = pca.transform(points)
	simplex = simplexpoints(c)
	updatedsphere = simplex+pca_x[0]
	updatepoints = pca.inverse_transform(updatedsphere)
	return updatepoints

def stereoproject(projpoint,projecteddsphere,center):
	plane = list(np.array(projpoint)-np.array(center))
	cofficient = sum(np.array(plane)*(np.array(center)))
	pts = []
	for i in range(0,len(projecteddsphere)):
		diff = projecteddsphere[i]-projpoint
		constterm =0;
		diffterm = 0;
		for x in range(0,len(plane)):
			constterm+= projpoint[x]*plane[x]
			diffterm+= diff[x]*plane[x]
		if(diffterm!=0):
			t = ((-1)*cofficient - constterm)/diffterm
			pnt = t*(projecteddsphere[i]-projpoint)+projpoint
			pts.append(pnt)
		else:
			pts.append(projecteddsphere[i])
	return PCAprojection(pts)

	
centroid = [0,0,0]
unitsphere = []
for x in newunitsphere:
	x = np.array(x)
	unitsphere.append(x/(math.dist(centroid,x)))
	
#xdata = [item[0] for item in unitsphere]
#ydata = [item[1] for item in unitsphere]
#zdata = [item[2] for item in unitsphere]
#ax.scatter3D(xdata, ydata, zdata,color = 'green')
#print("Rohit")		
#print(unitsphere)
#print("Rohit")		
	
#plt.show()
newfinalpoints = []

for y in unitsphere:
	plane = list(np.array(y)-np.array(centroid))
	cofficient = sum(np.array(plane)*(np.array(y)))
	dpts = [cofficient/x for x in plane]
	points = []
	i = 0
	points.append(y)
	for t in dpts[0:-1]:
		pnt = []
		for j in range(0,len(centroid)):
			if(i==j):
				pnt.append(t)
			else:
				pnt.append(0)
		points.append(pnt)
		i = i+1
	finalpoints = getsimplexsphere(points)
	newfinalpoints.extend(finalpoints.tolist())
projectionpoint = [random.randint(-1, 1) for x in range(len(newfinalpoints[0]))]
#print(projectionpoint)
projecteddsphere = []
for x in newfinalpoints:
		projecteddsphere.append(np.array(x)/math.dist([0,0,0],np.array(x)))
projpoint = np.array(projectionpoint)/math.dist([0,0,0],projectionpoint)
print(projpoint)


#print(len(projecteddsphere))
finalpoints = stereoproject(projpoint,projecteddsphere,centroid)

xdata = [item[0] for item in finalpoints]
ydata = [item[1] for item in finalpoints]

plt.scatter(xdata,ydata)
plt.savefig("projectionstereo.pdf", bbox_inches='tight')
plt.show()
ax = plt.axes(projection='3d')

xdata = [item[0] for item in projecteddsphere]
ydata = [item[1] for item in projecteddsphere]
zdata = [item[2] for item in projecteddsphere]
ax.scatter3D([projpoint[0]],[projpoint[1]],[projpoint[2]],color = 'green')
ax.scatter3D(xdata, ydata, zdata,color = 'blue')
plt.show()

'''
	c = len(centroid)-1
	pca = PCA(n_components=c).fit(points)
	pca_x = pca.transform(points)
	newpoints = []
	cnt = 0
	for x in range(0,c+1):
		coords = []
		for y in range(0,c+1):
			if(y==cnt):
				coords.append(1)
			else:
				coords.append(0)
		newpoints.append(coords)
		cnt = cnt+1
	dimlesspoints = PCAprojection(newpoints)
	summation = sum(dimlesspoints)
	dsphere = list(list(x) for x in dimlesspoints)
	for k in dimlesspoints:
		newcoord = (summation-k)/(len(dimlesspoints)-1)
		dsphere.append(list(newcoord))
	updatedsphere = []
	for pt in dsphere:
		updatedsphere.append(list(np.array(pt)+ pca_x[0]))
	dsphere.append(updatedsphere)
	updatedsphere = []
	print(dsphere)
	for lt in dsphere:
		lt = np.array(lt)
		print(lt)
		updatedsphere.append(lt/((math.dist([0,0],lt))*.01))
	#pca_x = np.append(pca_x,dsphere, axis=0)
	#pca_x.append([0,0])
	#print(updatedsphere)
	tr1 = [x12[0] for x12 in updatedsphere]
	tr2 = [x12[1] for x12 in updatedsphere]
	print(tr1,len(tr1))
	print(tr2,len(tr2))
	input()
	plt.plot(tr1,tr2)
	updatedsphere = np.array(updatedsphere)
	updatepoints = pca.inverse_transform(updatedsphere)
	newfinalpoints.append(list(updatepoints))
listr = list([list(tr1) for tr1 in tr] for tr in newfinalpoints)

print(len(listr))
zdata = [item[0] for item in listr]
xdata = [item[1] for item in listr]
ydata = [item[2] for item in listr]
ax.scatter3D(xdata, ydata, zdata,color = 'green')
plt.show()
print(listr)
'''
