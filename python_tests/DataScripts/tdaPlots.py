import numpy as np
from ripser import ripser
from persim import plot_diagrams
from sklearn import metrics
import argparse
import matplotlib as mpl
import pylab
from palettable.colorbrewer.qualitative import Set1_9
from palettable.tableau import Tableau_20
from persim import PersImage

##Arguments First
parser = argparse.ArgumentParser(description='Generate plots for given dataset')
parser.add_argument('--filename', type=str, help='Dataset filename', default='None')
parser.add_argument('--delimiter', type=str, help='Data delimiter', default=',')
parser.add_argument('--xlimit', help='Set a fixed limit for x-axis. (Default: use ceiling of max death)')
parser.add_argument('--nolegend', help='Do not print legends on plots.', action="store_true")
parser.add_argument('--inf', help='Value to replace infinite barcodes with (default use xlimit)')
args = parser.parse_args()

fig = mpl.pyplot.figure(figsize=mpl.pyplot.figaspect(.4))
#fig, (ax1, ax2) = mpl.pyplot.subplots(1,2)
#fig.set_size_inches(12,3.5)

if(args.filename == 'None'):
	data = np.random.random((100,2))
	fig.suptitle('Random data')
else:
	data = np.genfromtxt(args.filename, delimiter=args.delimiter)
	fig.suptitle(args.filename)

orig_data = data

## Generate the persistence intervals with ripser
diagrams = ripser(data,maxdim=2,thresh=5.0)['dgms']

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

ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

######BARCODES########

indexNaN = np.isinf(data)
if(args.inf):
	print "test"
	data[indexNaN] = float(args.inf)
else:
	data[indexNaN] = 1

# select our color map
colorPalette = Set1_9.mpl_colors
# set the markers for the scatter plotting of the persistence diagram.  
markers = ['o','+','<','d','x','>','1','2','3','4','^']

maxDeath = np.ceil(np.max(data[:,2]))
if args.xlimit is None :
    xLimit = maxDeath
else :
    xLimit = float(args.xlimit)

# did the user request an x-axis limit below the maximum death value found in the input?
if xLimit < maxDeath :
    print 'Requested xLimit (' + str(xLimit) + ') value is lower than max death value (' + str(maxDeath) + ')...aborting'
    sys.exit()

# are there more dimensions in the data then we have colors for?
if len(colorPalette) < np.max(data[:,0]) :
    print 'The current colormap has insufficient colors to represent all the dimensions in the data...aborting'
    sys.exit()

# build barcode plot

ax1.grid(True, which='both')
ax1.set_xlim(-.025,xLimit+.025)
ax1.set_xlabel('Time')
ax1.set_ylabel('Index')

for i in range(len(data)) :
    ax1.hlines(i+1, data[i,1], data[i,2], colors=colorPalette[int(data[i,0])])

# build the legend 
dimensions = []
for i in range(int(np.max(data[:,0]))+1) :
    dimensions.append(mpl.patches.Patch(color=colorPalette[i], label=r'$H_{}$'.format(i)))

#pylab.legend(handles=dimensions, loc='upper left', bbox_to_anchor=(.05,1))
if not args.nolegend :
    ax1.legend(handles=dimensions, loc='center right', bbox_to_anchor=(1,.65))

# build persistence diagram

axisMin = -0.025
axisMax = maxDeath+.025
ax2.grid(True, which='both')
ax2.set_xlim(axisMin, axisMax)
ax2.set_ylim(axisMin, axisMax)
ax2.set_xlabel('Birth')
ax2.set_ylabel('Death')
markerSize=15

# plot the data
for i in range(len(data)) :
    ax2.scatter(data[i,1], data[i,2], color=colorPalette[int(data[i,0])], marker=markers[int(data[i,0]) % len(markers)], s=markerSize)

# build the legend and build a 45 degree line
dimensions = []
for i in range(int(np.max(data[:,0]))+1) :
    dimensions.append(mpl.patches.Patch(color=colorPalette[i], label=r'$H_{}$'.format(i)))

if not args.nolegend :
    ax2.legend(handles=dimensions, loc='center right')

ax2.plot([axisMin, axisMax], [axisMin, axisMax], color='black', linestyle='-', linewidth=2)

pylab.savefig(args.filename + "_tdaPlot.pdf", bbox_inches='tight')
#pylab.show();

