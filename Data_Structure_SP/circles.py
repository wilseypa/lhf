import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib import cm
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

  
fig, ax = plt.subplots()

file1 = np.loadtxt("inputcentroid.txt",delimiter = " ")
x,y,r1 = file1[:,0],file1[:,1],file1[:,2]
patches = []
i = 0
col = ['cyan','orange','green','blue','red']
plt.scatter(x,y,s=5,c="red")
for x1, y1, r in zip(x, y, r1):
	if(r != 0):
		circle = Circle((x1, y1), r, edgecolor="violet",linewidth = 3,alpha = 0.4)
	else:
		ax.text(x1, y1, '({}, {})'.format(x1, y1))
		circle = Circle((x1, y1), radius = 0.01,edgecolor="red",linewidth = 1,alpha = 1)
	ax.add_patch(circle)
	i = i +1
ax.set_ylim([-1.2, 4])
ax.set_xlim([-2, 3])


l = []
file2 = np.loadtxt("simplices.txt",delimiter = " ")
for i in range(0,len(file2),3):
	if(file2[i:i+3,:].tolist() not in l):
		l.append(file2[i:i+3,:].tolist())
		t1 = plt.Polygon(file2[i:i+3,:], color="green",alpha=0.3)
		ax.add_patch(t1)
plt.show()


#file1 = np.loadtxt("inputcentroid.txt",delimiter = ",")
#x,y = file1[:,0],file1[:,1]
#plt.scatter(x,y)
#for i, j in zip(x, y):
#   plt.text(i, j, '({}, {})'.format(i, j))

#plt.show()