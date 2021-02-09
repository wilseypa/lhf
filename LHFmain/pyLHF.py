from Class import LHF
import tadasets
 

#Initialize the LHF Library  
pyLHF = LHF()


#Load data from file/generate data
data = tadasets.dsphere(n=100, d=3, r=1, noise=0.1)


#Load data into LHF to build array sizes (data types)
pyLHF.configData(data)


#Create argument dictionary and set any additional arguments
args = pyLHF.getDefaultArgs()


#Call a test function to ensure the library is properly loaded
pyLHF.testFunc(1)

#Call the runwrapper with default arguments and generated dataset
pyLHF.pyRunWrapper(args, data)

 
 
