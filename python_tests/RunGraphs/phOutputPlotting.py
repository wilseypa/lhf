
# Create barcode and persistent diagram plots of persistent homology outputs stored as a csv file of the form:

#     <dimension>, <birth>, <death>

# No assumption on the order of the rows in the csv file is assumed and the barcodes will be generated in the
# order that the data is stored in the input file.

# Support for plotting dimensions from 0-8 are provided.  The limiting factor on this is the palettable
# colors from Set1_9; simply using a palette with with a larger set of distinct colors will automatically
# extend this limit.

# The output files are written with input filename prefix appended with "-barcode.pdf" and "-pd.pdf".  The
# x-axis length can be set manually on the command line.  Program usage info is available with the command
# line option "--help" or "-h".

import sys
import os
import argparse
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import pylab
import numpy as np
from palettable.colorbrewer.qualitative import Set1_9
from palettable.tableau import Tableau_20

# here are some helper functions that we can use.
#--------------------------------------------------------------------------------

# define a function to display/save the pylab figures.
def displayGraph(fileName) :
    print("Creating graphics " + fileName)
    print("    ....writing pdf")
    pylab.savefig(fileName + ".pdf", bbox_inches='tight')
    pylab.clf()
#    pylab.show()
    return

####----------------------------------------------------------------------------------------------------
#### let us begin
####----------------------------------------------------------------------------------------------------

# select our color map
colorPalette = Set1_9.mpl_colors
# set the markers for the scatter plotting of the persistence diagram.  
markers = ['o','+','<','d','x','>','1','2','3','4','^']

# process the arguments on the command line
argparser = argparse.ArgumentParser(description='Create barcode and persistent diagram plots of persistent homology outputs from a csv file.')
argparser.add_argument('--xlimit', help='Set a fixed limit for x-axis. (Default: use ceiling of max death)')
argparser.add_argument('--nolegend', help='Do not print legends on plots.', action="store_true")
argparser.add_argument('--fileName', help='Name of the csv file.')
args = argparser.parse_args()

inFile = args.fileName
if inFile is None :
    print('Missing input filename....aborting')
    sys.exit()
    
rawData = np.loadtxt(inFile, dtype=np.float, delimiter=",", comments="#")

maxDeath = np.ceil(np.max(rawData[:,2]))
if args.xlimit is None :
    xLimit = maxDeath
else :
    xLimit = float(args.xlimit)

# did the user request an x-axis limit below the maximum death value found in the input?
if xLimit < maxDeath :
    print('Requested xLimit (' + str(xLimit) + ') value is lower than max death value (' + str(maxDeath) + ')...aborting')
    sys.exit()

# are there more dimensions in the data then we have colors for?
if len(colorPalette) < np.max(rawData[:,0]) :
    print('The current colormap has insufficient colors to represent all the dimensions in the data...aborting')
    sys.exit()

# build barcode plot

pylab.grid(True, which='both')
pylab.xlim(-.025,xLimit+.025)
pylab.xlabel('Time')
pylab.ylabel('Index')

for i in range(len(rawData)) :
    mpl.pyplot.hlines(i+1, rawData[i,1], rawData[i,2], colors=colorPalette[int(rawData[i,0])])

# build the legend 
dimensions = []
for i in range(int(np.max(rawData[:,0]))+1) :
    dimensions.append(mpl.patches.Patch(color=colorPalette[i], label=r'$H_{}$'.format(i)))

#pylab.legend(handles=dimensions, loc='upper left', bbox_to_anchor=(.05,1))
if not args.nolegend :
    pylab.legend(handles=dimensions, loc='center right', bbox_to_anchor=(1,.65))

displayGraph(inFile[0:len(inFile)-4] + '-barcode')

# build persistence diagram

axisMin = -0.025
axisMax = maxDeath+.025
pylab.grid(True, which='both')
pylab.xlim(axisMin, axisMax)
pylab.ylim(axisMin, axisMax)
pylab.xlabel('Birth')
pylab.ylabel('Death')
markerSize=60

# plot the data
for i in range(len(rawData)) :
    pylab.scatter(rawData[i,1], rawData[i,2], color=colorPalette[int(rawData[i,0])], marker=markers[int(rawData[i,0]) % len(markers)])

# build the legend and build a 45 degree line
dimensions = []
for i in range(int(np.max(rawData[:,0]))+1) :
    dimensions.append(mpl.patches.Patch(color=colorPalette[i], label=r'$H_{}$'.format(i)))

if not args.nolegend :
    pylab.legend(handles=dimensions, loc='center right')

pylab.plot([axisMin, axisMax], [axisMin, axisMax], color='black', linestyle='-', linewidth=2)

displayGraph(inFile[0:len(inFile)-4] + '-pd')
