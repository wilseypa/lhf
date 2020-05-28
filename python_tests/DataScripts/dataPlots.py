import numpy as np
from sklearn import metrics
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##Arguments First
parser = argparse.ArgumentParser(description='Generate plots for given dataset')
parser.add_argument('--filename', type=str, help='Dataset filename', default='None')
parser.add_argument('--delimiter', type=str, help='Data delimiter', default=',')
args = parser.parse_args()

fig = plt.figure(figsize=plt.figaspect(.3))
plt.axis('off')
#fig, (ax1, ax2) = mpl.pyplot.subplots(1,2)
#fig.set_size_inches(12,3.5)

if(args.filename == 'None'):
	data = np.random.random((100,2))
	plt.title('Random data')
else:
	data = np.genfromtxt(args.filename, delimiter=args.delimiter)
	plt.title(args.filename)

orig_data = data

####### Data plot #####

if(orig_data[0,:].size > 2):

    ax = fig.add_subplot(1,2,1,projection='3d')
    ax.scatter(orig_data[:,0],orig_data[:,1],orig_data[:,2])
    ax.set_xlabel('x0')
    ax.set_ylabel('x1')
    ax.set_zlabel('x2')

else:
    ax = fig.add_subplot(1,2,1)
    ax.scatter(orig_data[:,0],orig_data[:,1])

##### Distance matrix histogram #####

ax1 = fig.add_subplot(1,2,2)
dmat = metrics.pairwise_distances(orig_data);
dmat = np.reshape(dmat,(1,len(orig_data)**2))
dmat = np.sort(dmat)

print(dmat[0])

#ax = fig.add_subplot(1,4,4)
ax1.hist(dmat[0], bins=1000)#, cumulative = True)

plt.savefig(args.filename + "_dataplot.pdf", bbox_inches='tight')
#plt.show();

