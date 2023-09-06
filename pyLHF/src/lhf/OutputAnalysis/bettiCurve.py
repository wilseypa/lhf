#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
import matplotlib.patches as p
#import palettable
from matplotlib.cm import get_cmap
import math

def bettiCurve(perIntervals, fileName='', plotTitle = '', dimension=None, show=True, legend=True, hideTicks=False) :
    """
    Build betti curve for persistence intervals.

    Parameters
    ----------
    perIntervals : numpy.ndarray
        numpy [*,3] matrix; each row is a persistence interval written as <dimension, birth, death>.
    fileName : str
        name of the file for the persistence diagram pdf file.
    plotTitle : str
        title string for the persistence diagram (if '' no title generated, default: '')
    dimension : None or int
        dimension of betti curve to show; if none, the whole curve is plotted
    show : bool
        if true and the filename is not set the output is shown on a plot; otherwise the figure is returned to calling function
    legend : bool
        if true show the legend (H_d by color)
    hideTicks : bool
        if true hide the y-axis ticks (index of the PI)
    """
        
    
    bCurve = computeBettiCurve(perIntervals, dimension)
    
    ###   Plotting   ###
    
    #TODO: implement plotsize for figure layout; width > height (due to labels) 
    plotSize = 8.0
    plt.figure(figsize=(plotSize, plotSize/2))
    bCurve = np.array(bCurve)
    if dimension is None:
        plt.step(bCurve[:,0], bCurve[:,1], label='Betti Curve', where='post')
    else:
        plt.step(bCurve[:,0], bCurve[:,1], label=r'$H_{}$'.format(dimension) + ' Betti Curve', where='post')

    if legend:
        plt.legend()

    if plotTitle != '' :
        plt.set_title(plotTitle)

    if fileName != '' :
        plt.savefig(fileName, bbox_inches='tight')
    elif show:
        plt.show()
    else:
        return plt
    plt.close()

    return bCurve



    
def computeBettiCurve(perIntervals, dimension=None):
    """
    Compute the betti curve using dimensional persistence intervals. Setting a dimension only shows that
    dimension; if dimension is None (default) the Betti Curve / Euler Characteristic Curve is computed.
    Returns a list of lists representing the (0) filtration value and (1) betti number.

    Parameters
    ----------
    perIntervals : numpy.ndarray
        numpy [*,3] matrix; each row is a persistence interval written as <dimension, birth, death>.
    dimension : None or int
        dimension of betti curve to compute; if None, all dimensions are computed together (ECC)
        
    Returns
    ----------
    bettiCurve : list of tuples, [[epsilon, betti], ...]
        The Betti Curve or Dimensional Betti Curve for the persistence intervals.
    """
    bCurve = []
    bCurve.append([0.0, 0])

    for interval in perIntervals:
        #Each interval is [dim, birth, death]
        
        if(dimension != None):
            if(interval[0] != dimension or len(perIntervals) == 0):
                continue;
                
        #begin scanning to find the minimum value
        found = False
        ptr = 0
        lastResult = 0
        preSize = len(bCurve)

        #First loop to find birth time
        for z in range(len(bCurve)):
            if interval[1] == bCurve[z][0]: #Increment on found
                bCurve[z][1] += ((-1)**(interval[0]%2)*1)
                ptr = z + 1
                found = True
                lastResult = bCurve[z][1]
                break;

            if interval[1] < bCurve[z][0]: #Add if greater than
                bCurve.append([interval[1], lastResult + ((-1)**(interval[0]%2)*1) ])
                ptr = z
                found = True
                lastResult = bCurve[len(bCurve) - 1][1]
                break;

            lastResult = bCurve[z][1]

        if(not found):
            bCurve.append([interval[1], lastResult + ((-1)**(interval[0]%2)*1)])
            lastResult = bCurve[len(bCurve) - 1][1]
        
        ##Continue for death
        found = False
        for z in range(ptr, len(bCurve)):
            if interval[2] == bCurve[z][0]: #Decrement on found (deeath)
                bCurve[z][1] -= ((-1)**(interval[0]%2)*1)
                found = True
                ptr = z
                break;
                
            #Increment on value between
            elif interval[1] < bCurve[z][0] and interval[2] > bCurve[z][0]:
                bCurve[z][1] += ((-1)**(interval[0]%2)*1)
                ptr = z
                lastResult = bCurve[ptr][1]
            
            #Break on value above or end of list
            elif interval[2] > bCurve[z][0] or interval[2] < bCurve[z][0] or z == preSize - 1:
                break
                

        if(not found):
            bCurve.append([interval[2], lastResult - ((-1)**(interval[0]%2)*1)])

        #Sort the list (inefficient...)
        bCurve.sort(key=lambda x: x[0])

    return bCurve
