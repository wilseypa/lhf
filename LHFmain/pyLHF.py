from LHF import LHF
import sys
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
#pyLHF.runPH()

#pyLHF.runPH2()
#Call LHF with configured arguments and generated dataset
#store = bettiBoundaryTableEntry()
print ('Argument List:', str(sys.argv))
# print(type(sys.argv))
# print(type(sys.argv[1]))
# print(type(sys.argv[2]))
pyLHF.runPH3(sys.argv)
 
