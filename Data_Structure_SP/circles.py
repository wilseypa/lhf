
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
from random import randint
 
for epsilon in ["3.000000"]:#,"0.600000"]:#,"0.200000","0.300000","0.400000"]:
	for beta in ["0.100000"]:#,"1.000000","1.500000","3.000000"]:#,"0.050000","0.100000","0.200000","0.300000","0.400000","0.500000","0.600000","0.700000","0.800000","0.900000","1.000000","1.100000","1.200000","1.300000","1.400000","1.500000",]:
		filename = "simplicesE"+epsilon+"B"+ beta+".txt";
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		a = np.loadtxt(filename,delimiter = " ")
		#data = np.loadtxt("inputfile.txt",delimiter = ",")
		#ax.scatter(data[:,0],data[:,1],data[:,2],s=1,marker='o',c="black")

		ax.scatter(a[:,0],a[:,1],a[:,2],s=0.01,marker='o',c="black")

		a = np.reshape(a, (-1, 9))
		fc = ['#%06X' % randint(0, 0xFFFFFF) for i in range(a.shape[0])]

		poly3d = [[ a[i, j*3:j*3+3] for j in range(3)  ] for i in range(a.shape[0])]

		ax.add_collection3d(Poly3DCollection(poly3d, facecolors=fc, edgecolor = fc,linewidths=0.01,alpha = 1))
		plt.title(str(a.shape[0])+"-Simplices")
		#ax.set_xlim(-2,1.5)
		#ax.set_ylim(-1.5,2)
		#ax.set_zlim(0,3)
		ax.view_init(elev=-107, azim=25)
		plt.axis('off')
		plt.savefig("plot"+filename+".pdf",bbox_inches='tight',pad_inches = 0)
		plt.show()
		plt.close()
'''

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
		#ax.text(x1, y1, '({}, {})'.format(x1, y1))
		circle = Circle((x1, y1), radius = 0.01,edgecolor="red",linewidth = 1,alpha = 1)
	ax.add_patch(circle)
	i = i +1

l = []
file2 = np.loadtxt("simplices.txt",delimiter = " ")
print(file2[0:5])
for i in range(0,len(file2),3):
	if(file2[i:i+3,:].tolist() not in l):
		l.append(file2[i:i+3,:].tolist())
		t1 = plt.Polygon(file2[i:i+3,:], color="green",alpha=0.3)
		ax.add_patch(t1)
plt.savefig("MeshGenerationSP2.pdf", bbox_inches='tight')
plt.show()


#file1 = np.loadtxt("inputcentroid.txt",delimiter = ",")
#x,y = file1[:,0],file1[:,1]
#plt.scatter(x,y)
#for i, j in zip(x, y):
#   plt.text(i, j, '({}, {})'.format(i, j))

#plt.show()
'''
