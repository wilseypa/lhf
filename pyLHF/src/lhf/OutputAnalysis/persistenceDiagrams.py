#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
#import palettable
from matplotlib.cm import get_cmap
import math

#### This is the basic code to plot persistence diagrams.  We could go crazy with programming options, but I don't recomment it.
#### The one thing we might want to add is the ability to add histograms on the top and left of the quadrant to show point
#### populations in a specific homology group (I suppose it could be done for multiple homology groups, but I don't think that
#### would be wise [or worth the effort])

def persistenceDiagram(perIntervals, fileName='', plotTitle = '', maxDist=-1.0, projPlts='', dimRange=[], numBins=100, show=True) :
    """
    Build Persistence Diagrams for persistence intervals.

    Build a matplotlib persistence diagram (scatter plot) for the input persistence intervals (perIntervals).  The plotTitle
    string will be placed on the plot if it exists.  The plot will be displayed on the computer screen unless a fileName parameter
    is defined (in which case, the plot will be written to the file).  The suffix on the fileName string will control the format
    (pdf, png, svg, ...) of the written file.

    A legend will be generated and placed in the lower right corner of the plot.  The legend will use homology group numbering
    (H0, H1, ...) for the dimensions and the number of points in each homology group will be written with parenthesis in the
    legend. 

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
    projPlts : str 
        generate x and y axis projection plots (from a histogram); possible plot types: '' (none), 'counts', 'density'
    dim : [ int ]
        list of dimensions to show in projection plots; if empty, all dimensions shown
    numBins : int 
        number of bins in the histogram plots (default=100)
    show : bool
        if true and the filename is not set the output is shown on a plot; otherwise the figure is returned to calling function
    """

    # verifty legal values for numPlts; if not, turn off
    if not projPlts in ['counts', 'density'] :
        projPlts = ''

    # if no dimensions requested, plot all
    if not dimRange :
        dimRange = range(int(np.max(perIntervals[:,0])+1))

    histAsDensity = False
    if projPlts == 'density' :
        histAsDensity = True

    infInPIs = False
    maxDeath = np.max(perIntervals[perIntervals[:,2] != np.inf][:,2])

    # if inf exists in the death of a PI, expand the max axis by .01 and set the inf values to this
    if not np.isfinite(perIntervals[:,2]).all() :
        infInPIs = True
        maxDeath = maxDeath + (.01 * maxDeath)
        np.nan_to_num(perIntervals, posinf=maxDeath, copy=False)
        
    if maxDist == -1.0 or maxDeath > maxDist:
        maxDist = maxDeath

    # s: list of scatter plots for input to the plt.legend function
    s = []
    # lbls: hold dimension labels (H_0, H_1, ...)
    lbls = []

    #colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors
    colors = get_cmap('tab20').colors
    markers = ['.', '1', '2', '3', '4', '+', 'x']
    colorIndx = 0
    markerIndx = 0

    fontSize = 14
    plt.rcParams.update({'font.size': fontSize})

    # definitions for the axes
    plotSize = 8.0
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    
    if projPlts != '' :
        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]
    else :
        rect_scatter = [left, bottom, .8, .8]

    plt.figure(figsize=(plotSize, plotSize))

    ## also shamelessly copied....
    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    if projPlts != '' :
        ax_histx = plt.axes(rect_histx)
        ax_histx.tick_params(direction='in', labelbottom=False)
        ax_histy = plt.axes(rect_histy)
        ax_histy.tick_params(direction='in', labelleft=False)
        ax_histx.set_xlim(0, maxDist+(.02*maxDist))
        ax_histy.set_ylim(0, maxDist+(.02*maxDist))

    ax_scatter.axis('equal')

    if plotTitle != '' :
        ax_scatter.set_title(plotTitle)

    ax_scatter.set_xlim(0, maxDist+(.02*maxDist))
    ax_scatter.set_ylim(0, maxDist+(.02*maxDist))
    
    ax_scatter.set_ylabel('Death')
    ax_scatter.set_xlabel('Birth')

    # plot black diagonal line
    ax_scatter.plot(np.arange(maxDist+1), color='black', linewidth=.6)

    # plot horizontal line for infinity if value present in PIs
    if infInPIs :
        ax_scatter.axhline(y = maxDeath, color = 'r', linestyle = 'dashed', label='inf')
        #ax_scatter.annotate(r'$\infty$', xy=(0,maxDeath), xytext=(0-(.06*maxDeath),maxDeath-(.02*maxDeath)), fontsize=26)

    # my solution for building the legends may be more complicated than necessary, but i built it and it works, so i'm going to
    # leave it as is

    # plot each dimension in order (so the higher dimensions will be on top)
    #for i in range (math.floor(np.min(perIntervals[:,0])), math.ceil(np.max(perIntervals[:,0])+1)) :
    for i in dimRange :
        # if no persistence intervals in this dimension skip to the next
        if len(perIntervals[:, [0]][perIntervals[:, 0] == i]) == 0 : continue
        x = perIntervals[:, [1]][perIntervals[:, 0] == i] # birth at dimension i
        y = perIntervals[:, [2]][perIntervals[:, 0] == i] # death at dimension i
        s.append(ax_scatter.scatter(x, y, c = np.atleast_2d(colors[colorIndx]), marker=markers[markerIndx], s=60))
        lbls.append(r'${H_{%s}}$' % '{:,}'.format(i) + ' (' + str(len(perIntervals[:, [1]][perIntervals[:, 0] == i])) + ')' )
        if (projPlts != '') :
            ax_histx.hist(x, bins=numBins, range=(0, maxDist), color=np.atleast_2d(colors[colorIndx]), density=histAsDensity)
            ax_histy.hist(y, bins=numBins, range=(0, maxDist), color=np.atleast_2d(colors[colorIndx]), density=histAsDensity, orientation='horizontal')
        colorIndx += 1
        if colorIndx == 20 :
            colorIndx = 0
            markerIndx += 1

    ax_scatter.legend(reversed(s), reversed(lbls), loc='lower right')

    if fileName != '' :
        plt.savefig(fileName, bbox_inches='tight')
    elif show:
        plt.show()
    else:
        return plt
    plt.close()

    return
