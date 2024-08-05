#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
import math


def show_values(pc, fmt="%.4f", **kw):
    """
    Adds value labels to the heatmap for the dissimilarity matrix


    Parameters
    ----------
    pc : numpy.ndarray
        numpy square matrix with dissimilarity values
    fmt : str
        format of decimal values on the plot
    """
    
    pc.update_scalarmappable()
    ax = pc.axes
    for p, color, value in zip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
        x, y = p.vertices[:-2, :].mean(0)
        if np.all(color[:3] > 0.5):
            color = (0.0, 0.0, 0.0)
        else:
            color = (1.0, 1.0, 1.0)
        ax.text(x, y, fmt % value, ha="center", va="center", color=color, **kw)
        
        

def heatmap(disMatrix, fileName='', plotTitle = '', columns=[], show=True, datalabel=True) :
    """
    Build heatmap diagram for dissimilarity matrix


    Parameters
    ----------
    disMatrix : numpy.ndarray
        numpy square matrix with dissimilarity values
    fileName : str
        name of the file for the persistence diagram pdf file.
    plotTitle : str
        title string for the heatmap (if '' no title generated, default: '')
    show : bool
        if true and the filename is not set the output is shown on a plot; otherwise the figure is returned to calling function
    """

    #TODO: implement plotsize for figure layout; width > height (due to labels) 
    #plotSize = 8.0
    #plt.figure(figsize=(plotSize, plotSize))
    
    fig, ax = plt.subplots()
    c = ax.pcolor(disMatrix, edgecolors='k', linestyle='dashed', linewidths=0.2, cmap='Blues')#, vmin=0.0, vmax=.04)

    ax.set_yticks(np.arange(disMatrix.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(disMatrix.shape[1]) + 0.5, minor=False)
    ax.set_xticklabels(columns, minor=False)
    ax.set_yticklabels(columns, minor=False)
    #ax.xaxis.tick_top()

    if plotTitle != '' :
        ax_scatter.set_title(plotTitle)

    plt.xlim((0,disMatrix.shape[1]))

    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    plt.gca().invert_xaxis()

	# Datalabel outputs values on each square in the dissimilarity matrix
    if(datalabel):
        show_values(c)
        
    if fileName != '' :
        plt.savefig(fileName, bbox_inches='tight')
    elif show:
        plt.show()
    else:
        return plt
    plt.close()

    return
