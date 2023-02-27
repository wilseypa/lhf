import sys
sys.path.append('../')
import tadasets
import numpy as np

from LHF import LHF
#from LHF import pipePacket
#from LHF import bettiBoundaryTableEntry 

#Load data from file/generate data
data = tadasets.dsphere(n=100, d=3, r=1, noise=0.1)
np.savetxt("tempData.csv", data, fmt='%.4f', delimiter=',')



#Initialize the LHF Library  
pyLHF = LHF.LHF()

#Set debug mode to true, configure other arguments (optional)
pyLHF.args["debug"] = "1"
pyLHF.args["inputFile"] = "tempData.csv"
pyLHF.args["dimensions"] = "3"
pyLHF.args["complexType"] = "simplexArrayList"
pyLHF.args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence"
pyLHF.args["epsilon"] = "5.0"

#Call a test function to ensure the library is properly loaded
pyLHF.testFunc(1, "HELLO")


pyLHF.runPH(data)
 
