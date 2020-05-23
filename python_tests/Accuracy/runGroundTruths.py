import numpy as np
from ripser import ripser
import argparse
import sys

##Arguments First
parser = argparse.ArgumentParser(description='Generate k-means++ reduced persistence intervals for given dataset')
parser.add_argument('--filename','-f', type=str, help='Dataset filename', default='None')
parser.add_argument('--dim','-d', type=int, help='Homology dimension', default=1)
parser.add_argument('--epsilon','-e',type=float, help='Epsilon threshold', default=2.0)
args = parser.parse_args()

if(args.filename == 'None'):
	print "No filename specified... returning"
	sys.exit(0)
else:
	data = np.genfromtxt(args.filename, delimiter=',')

print args, len(data)

## Generate the persistence intervals with ripser
diagrams = ripser(data,maxdim=args.dim,thresh=args.epsilon)['dgms']

## Do some sorting and formatting of the persistence intervals
persArray = []

i = 0;
for dim in diagrams:
	tempDim = [];
	for interval in dim:
		interval = np.insert(interval, 0, i)
		tempDim.append(interval)
		
	tempDim.sort(key=lambda x: x[1])
	i+=1
	
	persArray = persArray + tempDim
	
data = np.asarray(persArray)

indexNaN = np.isinf(data)
	
data[indexNaN] = float(args.epsilon)

np.savetxt(args.filename + "_d" + str(args.dim) + "_PH", data, delimiter=',')
