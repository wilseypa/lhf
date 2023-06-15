import sys
sys.path.append('../')
import tadasets
import persim
import numpy as np
import matplotlib.pyplot as plt

from LHF import LHF
from LHF.OutputAnalysis import persistenceDiagram, heatmap, barcodeDiagram, bettiCurve

from LHF.DataGeneration import dataGen as dg
from LHF.DataGeneration import objGen as og

pis = []



#def fish(param):
#	a= 1
#	return (param[0]**2+param[1]**2)**2 - a*param[0]*(param[0]**2 - param[1]**2) 

def plot(pointCloud) :

    plt.title("Point Cloud, Total Points: %s" % pointCloud.shape[0])
    plt.scatter(pointCloud[:,0],pointCloud[:,1], marker=".", s=1)
    plt.axis('square')
    plt.show()

    return


for i in range(0, 1):


    #Initialize the LHF Library  
    pyLHF = LHF.LHF()

    #Set debug mode to true, configure other arguments (optional)
    pyLHF.args["debug"] = "1"
    pyLHF.args["mode"] = "reduced"
    pyLHF.args["clusters"] = 200
    pyLHF.args["dimensions"] = 2
    
    #Generate a square of points and make cuts using new dataGen library
    #square = dg.genFilledCube(dim=2)
    
    #d = dg.buildObj('(_x[0]**2) + _x[1]**2 <= 0.8', square)
    #d = dg.buildObj('\
    #_x[0] < -0.25 or \
    #(_x[0]+0.25)**2 + (_x[1]-0.25)**2 > 0.6**2\
    #', d)
    
    
    #Fish example
    d = og.butterflyCurve(3000, 2)
    np.savetxt('./pointCloud.csv',d,delimiter=',')
    plot(d)
    
    #Run PH and get the full bettiTable, pipePacket object
    boundpis, ppkt, elapsed = pyLHF.runPH(d)
       
    ##Extract the raw pis from the bettiTable (i.e. remove last column, boundary generators)
    pis = np.array([[z[0],z[1],z[2]] for z in boundpis])
        


plt = bettiCurve(pis, dimension=1, show=True)


