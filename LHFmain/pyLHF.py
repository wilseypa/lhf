from LHF import LHF
#from LHF import bettiBoundaryTableEntry 
import tadasets

#Load data from file/generate data
data = tadasets.dsphere(n=100, d=3, r=1, noise=0.1)

#Initialize the LHF Library  
pyLHF = LHF(data)

#Set debug mode to true, configure other arguments (optional)
pyLHF.args["debug"] = "1"
pyLHF.args["dimensions"] = "3"

#Call a test function to ensure the library is properly loaded
pyLHF.testFunc(1, "HELLO")

#Call LHF with configured arguments and generated dataset
pyLHF.runPH()


#Call LHF with configured arguments and generated dataset
#store = bettiBoundaryTableEntry()
pyLHF.runPH2()
 
