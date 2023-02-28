import sys
sys.path.append('../')
import tadasets
import persim
import numpy as np

from LHF import LHF
from LHF.OutputAnalysis import persistenceDiagram, heatmap, barcodeDiagram, bettiCurve

#Load data from file/generate data
data = tadasets.dsphere(n=100, d=3, r=1, noise=0.1)
np.savetxt("tempData.csv", data, fmt='%.4f', delimiter=',')



#Initialize the LHF Library  
pyLHF = LHF.LHF()

#Set debug mode to true, configure other arguments (optional)
pyLHF.args["debug"] = "0"
pyLHF.args["epsilon"]= 1.0

pis = pyLHF.runPH(data)

 
print(pis)

plt = persistenceDiagram(pis, show=True)


