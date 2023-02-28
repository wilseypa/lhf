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
pyLHF.args["debug"] = "0"


pis = pyLHF.runPH(data)
 
print(pis)
