import sys
sys.path.append('../')
import tadasets
import persim
import numpy as np
import matplotlib.pyplot as plt


import lhf

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
    pyLHF = lhf.createPipeline()

    #Set debug mode to true, configure other arguments (optional)
    pyLHF.config["debug"] = "1"
    pyLHF.config["mode"] = "reduced"
    pyLHF.config["clusters"] = 200
    pyLHF.config["dimensions"] = 2
    
    #Generate a square of points and make cuts using new dataGen library
    cube = lhf.DataGeneration.genFilledCube(dim=3)
    
    #Cut into a filled sphere
    #d = dg.buildObj('_x[0]**2 + _x[1]**2 + _x[2]**2 <= 1', cube)
    
    #Spiral Example
    #d = dg.buildObj('_x[0]**2 + _x[1]**2 + _x[2]**2 <= 1', cube)
    #d = dg.buildObj('(_x[0]-0.7*math.cos(_x[2]/.25))**2 + (_x[1]-0.7*math.sin(_x[2]/.25))**2 <= 0.9', d)
    #d = dg.buildObj('(_x[0]**2)+(_x[1]**2) >= 0.1', d)
    
    #Square Spiral Example
    #d = dg.buildObj('(_x[0]-0.7*math.cos(_x[2]/.25))**2 + (_x[1]-0.7*math.sin(_x[2]/.25))**2 <= 0.9', cube)
    #d = dg.buildObj('(_x[0]**2)+(_x[1]**2) >= 0.1', d)
    
    #Half bowl with sphere
    #d = dg.buildObj('(_x[2] < 0 and ( _x[0]**2 + _x[1]**2 + (_x[2]+.25)**2 <= 1 and _x[0]**2 + _x[1]**2 + (_x[2]+.25)**2 >= 0.55)) or \
    #                (_x[0]**2 + _x[1]**2 + (_x[2]-.5)**2 <=0.35 and _x[0]**2 + _x[1]**2 + (_x[2]-.5)**2 >=0.15)', cube)
    
    
    #Mushroom/stamp
    #d = dg.buildObj('(_x[2] <= 0.05 and ( _x[0]**2 + _x[1]**2 + (_x[2]+.25)**2 <= 1 and (_x[0]**2 + _x[1]**2 + (_x[2]+.25)**2 >= 0.55 or (abs(_x[2]) <=0.05 and _x[0]**2 + _x[1]**2 <= 1)))) or \
    #                (_x[0]**2 + _x[1]**2 + (_x[2]-.5)**2 <=0.35 and _x[0]**2 + _x[1]**2 + (_x[2]-.5)**2 >=0.15)', cube)
    
    #Spiral tube
    d = lhf.DataGeneration.buildObj('(_x[0]-.5*math.cos(_x[2]/.2))**2 + (_x[1]-.5*math.sin(_x[2]/.15))**2 <= 1.0 and \
                    (_x[0]-.5*math.cos(_x[2]/.2))**2 + (_x[1]-.5*math.sin(_x[2]/.15))**2 >= 0.1', cube)
    
    #Cut a torus
    #d = dg.buildObj('(_x[0]**2 + _x[1]**2)**.5 + _x[2]**2 <= 0.5', cube)
    
    #Cut a half hole through the center of the sphere
    #d = dg.buildObj('_x[0] > 0 or _x[1]**2 + _x[2]**2 >= 0.05', d)
    
    #Cut a void in the sphere
    #d = dg.buildObj('_x[0]**2 +_x[1]**2 + _x[2]**2 >=0.7 or (abs(_x[0]) < 0.05)', d)
    
    
    #Attempt to cut a spiral
    #d = dg.buildObj('(_x[0]-0.5*math.cos(_x[2]/.25))**2 + (_x[1]-0.5*math.sin(_x[2]/.25))**2 <= 0.9 and \
    #                (_x[0]-0.5*math.cos(_x[2]/.05))**2 + (_x[1]-0.5*math.sin(_x[2]/.05))**2 >= 0.6', d)
    
    
    
    #d = dg.buildObj('(_x[0]**2) + _x[1]**2 <= 0.8', square)
    #d = dg.buildObj('\
    #_x[0] < -0.25 or \
    #(_x[0]+0.25)**2 + (_x[1]-0.25)**2 > 0.6**2\
    #', d)
    
    
    #Fish example
    #d = og.butterflyCurve(3000, 2)
    np.savetxt('./pointCloud.csv',d,delimiter=',')
    #plot(d)
    
    exit(0)
    
    #Run PH and get the full bettiTable, pipePacket object
    boundpis, ppkt, elapsed = pyLHF.runPH(d)
       
    ##Extract the raw pis from the bettiTable (i.e. remove last column, boundary generators)
    pis = np.array([[z[0],z[1],z[2]] for z in boundpis])
        


plt = bettiCurve(pis, dimension=1, show=True)


