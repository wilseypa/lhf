import sys
sys.path.append('../')
import tadasets
import persim
import numpy as np

from LHF import LHF
from LHF.OutputAnalysis import persistenceDiagram, heatmap, barcodeDiagram, bettiCurve




#Initialize the LHF Library  
pyLHF = LHF.LHF()

#Set debug mode to true, configure other arguments (optional)
pyLHF.args["debug"] = "0"
pyLHF.args["epsilon"]= 1.0

pis = []

for i in range(0, 10):
	
	#Load data from file/generate data
	data = tadasets.dsphere(n=100+(i*15), d=3, r=1, noise=0.1)

	pis = pyLHF.runPH(data)

	print(len(pis))



plt = persistenceDiagram(pis, show=True)


