#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
import matplotlib.patches as p
#import palettable
from matplotlib.cm import get_cmap
from matplotlib import gridspec
from scipy.spatial import distance_matrix
import math
from barcodeDiagrams import barcodeDiagram

    

def triangulatedPlot(origData, ax=None, eps=0.0, cmap='tab10'):
    """
    Build triangulated plot at a specified epsilon value
    
    Parameters
    ----------
    origData : numpy.ndarray
        numpy [*,*] matrix; each row is a sample of dimension d; only dimensions d0 and d1 are used for triangulation plot
    ax : matplotlib.axes
        if set the plot is drawn to the passed axis
    eps : float
        epsilon value to draw the triangulated point cloud at
    cmap : str
	    color map used for coloring different dimensional components between barcode plot and triangulated plot
    """
    
    if(ax == None):
        fig, ax = plt.subplots()
        
    #Get edges <= eps
    edges = []
    dmat = distance_matrix(origData, origData)
    
    for i in range(len(dmat)):
        for j in range(i+1, len(dmat)):
            if dmat[i,j] < eps:
                edges.append([i, j])
                
    colors = get_cmap(cmap).colors
                
    #Get triangles <= eps
    tris = []
    
    for e1 in range(len(edges)):
        for e2 in range(e1+1, len(edges)):
            for e3 in range(e2+1, len(edges)):
                nbrs = edges[e1] + edges[e2] + edges[e3]
                if(len(set(nbrs)) == 3):
                    tris.append([*set(nbrs)])
    
    #Plot edges and triangles
    ax.scatter(origData[:,0], origData[:,1], s=3, c='black')
    for edge in edges:
        pt1, pt2 = [origData[pt] for pt in [n for n in edge]]
        line = plt.Polygon([pt1, pt2], closed=None, fill=None, edgecolor=colors[0])
        plt.gca().add_patch(line)

    for tri in tris:
        pt1, pt2, pt3 = [origData[pt] for pt in [n for n in tri]]
        line = plt.Polygon([pt1,pt2,pt3], closed=False, color=colors[1], alpha=0.01, fill=True, edgecolor=None)
        plt.gca().add_patch(line)

    ax.set_yticks([])
    ax.set_xticks([])
    ax.axis('off')
    
    ax.set_title('$\epsilon$ = ' + str(eps))
    
    return



def filtrationImage(origData, perIntervals, fileName='', plotTitle = '', maxDist=-1.0, show=True, cmap='tab10') :
    """
    Build Filtration Image (triangulated point cloud + barcode plot) for input data set and persistence intervals

    Parameters
    ----------
    origData : numpy.ndarray
        numpy [*,*] matrix; each row is a sample of dimension d; only dimensions d0 and d1 are used for triangulation plot
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
    cmap : str
	    color map used for coloring different dimensional components between barcode plot and triangulated plot
    """
    
    fig = plt.figure(figsize=(10.0, 6.0))
    gs = gridspec.GridSpec(2, 1, figure=fig)
    
    #Build grid spec for triangulated images
    gs2 = gs[0].subgridspec(1,4)
    
    ax0 = fig.add_subplot(gs2[0,0])
    triangulatedPlot(origData, ax=ax0, eps=0.0, cmap=cmap)
    
    ax1 = fig.add_subplot(gs2[0,1])
    triangulatedPlot(origData, ax=ax1, eps=0.5, cmap=cmap)
    
    ax2 = fig.add_subplot(gs2[0,2])
    triangulatedPlot(origData, ax=ax2, eps=1.0, cmap=cmap)
    
    ax3 = fig.add_subplot(gs2[0,3])
    triangulatedPlot(origData, ax=ax3, eps=1.5, cmap=cmap)

    ax4 = fig.add_subplot(gs[1])
    barcodeDiagram(perIntervals, show=False, ax=ax4, cmap=cmap)
    
    ax4.set_xlabel('$\epsilon$')
    ax4.set_yticks([])

    plt.tight_layout()
    
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


if __name__ == "__main__":
    
    import argparse
    from sklearn.datasets import make_circles
    from ripser import ripser
    
    ##Arguments First
    parser = argparse.ArgumentParser(description='Generate Filtration Image')
    parser.add_argument('--filename', '-f', type=str, help='Input filename (default: makeCircles)', default=None)
    parser.add_argument('--output', '-o', type=str, help='Output image filename (default: show plot)', default='')
    parser.add_argument('--noPoints', '-n', type=int, help='Number of points if generating makeCircles', default=40)
    parser.add_argument('--plotTitle', '-p', type=str, help='Title for the plot (default: \'\')', default='')
    parser.add_argument('--cmap', '-c', type=str, help='Matplotlib cmap to use for coloring (default: tab10)', default='tab10')
    args = parser.parse_args()
    
    if(args.filename):
        data = np.loadtxt(args.filename, delimiter=',')
    else:
        data = make_circles(args.noPoints, noise=0.15)[0]
        np.savetxt('makeCircles_' + str(args.noPoints) + '.csv', data, delimiter=',')
    
    diagrams = ripser(data, maxdim=len(data[0])-1)['dgms']
    fmtDiagrams = []
    
    for dim in range(len(diagrams)):
        for i in diagrams[dim]:
            fmtDiagrams.append([dim, i[0], i[1]])
    
    fmtDiagrams = np.array(fmtDiagrams)
    
    filtrationImage(data, fmtDiagrams, args.output, args.plotTitle, show=True, cmap=args.cmap)

    
    
