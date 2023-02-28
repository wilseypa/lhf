#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
import matplotlib.patches as p
#import palettable
from matplotlib.cm import get_cmap
import math

def barcodeDiagram(perIntervals, fileName='', plotTitle = '', maxDist=-1.0, show=True, legend=True, hideTicks=False, ax=None, cmap='tab10') :
    """
    Build Barcode Diagrams for persistence intervals.

    Parameters
    ----------
    perIntervals : numpy.ndarray
        numpy [*,3] matrix; each row is a persistence interval written as <dimension, birth, death>.
    fileName : str
        name of the file for the persistence diagram pdf file.
    plotTitle : str
        title string for the persistence diagram (if '' no title generated, default: '')
    maxDist : float
        x/y axis upper bound of the persistence diagram, defaults to -1 (use the max birth/death value). This value will be replaced with the range of birth/death limits if larger.
    show : bool
        if true and the filename is not set the output is shown on a plot; otherwise the figure is returned to calling function
    legend : bool
        if true show the legend (H_d by color)
    hideTicks : bool
        if true hide the y-axis ticks (index of the PI)
    ax : matplotlib.axes
        if set the plot is drawn to the passed axis
    cmap : str
	    color map used for coloring different dimensional components between barcode plot and triangulated plot
    """

    #TODO: implement plotsize for figure layout; width > height (due to labels) 
    #plotSize = 8.0
    #plt.figure(figsize=(plotSize, plotSize))
    if(ax == None):
        fig, ax = plt.subplots()

    infInPIs = False
    maxDeath = np.max(perIntervals[perIntervals[:,2] != np.inf][:,2])

    # if inf exists in the death of a PI, expand the max axis by .01 and set the inf values to this
    if not np.isfinite(perIntervals[:,2]).all() :
        infInPIs = True
        maxDeath = maxDeath + (.01 * maxDeath)
        np.nan_to_num(perIntervals, posinf=maxDeath, copy=False)
        
    if maxDist == -1.0 or maxDeath > maxDist:
        maxDist = maxDeath
    
    #Perform sorting of the intervals for visual appeal
    perIntervals = perIntervals[perIntervals[:,2].argsort()]
    perIntervals = perIntervals[perIntervals[:,1].argsort(kind='mergesort')]
    perIntervals = perIntervals[perIntervals[:,0].argsort(kind='mergesort')]
    
    colors = get_cmap(cmap).colors

    xLimit = float(maxDeath)

    # build barcode plot
    ax.grid(True, which='both')
    ax.set_xlim(-.025,xLimit+.025)
    
    if hideTicks:
        ax.set_yticks([])

    for i in range(len(perIntervals)) :
        ax.hlines(i+1, perIntervals[i,1], perIntervals[i,2], colors=colors[int(perIntervals[i,0])])

    # build the legend 
    dimensions = []
    for i in range(int(np.max(perIntervals[:,0]))+1) :
        dimensions.append(p.Patch(color=colors[i], label=r'$H_{}$'.format(i)))

    if legend and ax:
        ax.legend(handles=dimensions, loc='lower right')
    elif legend:
        ax.legend(handles=dimensions, loc='lower right', bbox_to_anchor=(1,.65))
		

    if plotTitle != '' :
        ax.set_title(plotTitle)

    if fileName != '' :
        plt.savefig(fileName, bbox_inches='tight')
    elif show:
        plt.show()
    else:
        return plt
    plt.close()

    return
