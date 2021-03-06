import sys
import numpy as np
import argparse
import csv
import subprocess
import tadasets
import time
import os
import string
np.set_printoptions(threshold=sys.maxsize) # for debugging
from sklearn.datasets import make_blobs



''' MAIN LOOP HERE '''

''' Get CMD Arguments: [filename] [epsilon] [maxDim] [partitioner] [upscale 0/1] [centroids #] [reps]'''
parser = argparse.ArgumentParser(description='Run iterative testing for TDA tools')
parser.add_argument('--epsilon','-e',type=float, help='Max epsilon to compute PH up to', default=5)
parser.add_argument('--dim', '-d', type=int, help='Maximum homology dimension to compute (Hx)', default=1)
parser.add_argument('--datadim', '-nd', type=int, help='Dimension of data to generate', default=2)
parser.add_argument('--npoints', '-n', type=int, help='Number of points to generate', default=100)
parser.add_argument('--centroids', '-k', type=int, help='Number of centroids to compute', default=50)
parser.add_argument('--np','-np', type=int, help='Number of processors to use', default=1)
parser.add_argument('--reps','-r',type=int,help='Number of experiment repititions', default = 1)
parser.add_argument('--mode','-m',type=str,help='Mode to run LHF in', default='mpi')
args = parser.parse_args()

epsilon = args.epsilon
maxDim = args.dim
dataDim = args.datadim
npoints = args.npoints
centroids = args.centroids
reps = args.reps
mode = args.mode
noproc = args.np
timings = {}

print "Running " + mode + " experiments with:\nEpsilon: ", epsilon, "\nMaxDim", maxDim, "\ndataDim: ",  dataDim, "\nnPoints: ", npoints, "\nCentroids: ", centroids, "\nNProcs", noproc, "\n"

''' Generate input data '''
originalData = tadasets.dsphere(n=npoints, d=dataDim-1, r=1, noise=0.1)

for i in range(0, int(reps)):

	''' Create folder for output '''
	outDir = "perf_dsphere_"+str(centroids)+"_e"+str(epsilon) + "_d"+str(maxDim) + "_np" + str(noproc) +"/" + str(i) + "/"
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	if not os.path.isfile(os.getcwd() + "/aggResults.csv"):
		outfile = file(os.getcwd() + "/aggResults.csv", 'a')
		outfile.write('OutputPath,No.Procs,Vectors,Dimensions,maxDim,Epsilon,Centroids,PH_Time(s),BettiCount\n')
		outfile.close()      
        
        
	''' Output all data to files '''
	np.savetxt(outDir + "originalData.csv", originalData, delimiter=',')

	#Check if the last line was returned, keeping our data clean
	isRet = True;
	with file(os.getcwd() + "/aggResults.csv", 'rb+') as f:
		f.seek(-1,2)
		a = f.read()
		if(a is not "\n"):
			isRet = False;
	
	outfile = file(os.getcwd() + "/aggResults.csv", 'a')
	if(not isRet):
		outfile.write('')
	outfile.write(outDir[:-1] + "," + str(noproc) + "," + str(npoints) + "," + str(dataDim) + "," + str(maxDim) + "," + str(epsilon) + "," + str(centroids) + ",")
	outfile.close()

	
	#try:
	start = time.time()
	print "mpiexec","-np",str(noproc),"./LHF", "-m",mode,"-d", str(maxDim), "-e", str(epsilon), "-k",str(centroids), "-s", "1", "-i", os.getcwd()+"/"+outDir + "originalData.csv"
	proc = subprocess.check_output(["mpiexec","-np",str(noproc),"./LHF", "-m",mode,"-d", str(maxDim), "-e", str(epsilon), "-k",str(centroids), "-s", "1", "-i", os.getcwd()+"/"+outDir + "originalData.csv"])  


	end = time.time()
	pers_time = (end - start)
	outdata = np.genfromtxt("./output.csv.csv")
	
	outfile = file(os.getcwd() + "/" + outDir +"LHF_Output.csv", 'w')
			
	outfile.write(outdata)
			
	outfile.close()
			
	outfile = file(os.getcwd() + "/aggResults.csv", 'a')
	outfile.write(str(pers_time) + "," + str(outdata.size - 1) + ",")
	outfile.close()
   #except:
		#print "LHF processing failed"    


	outfile = file(os.getcwd() + "/aggResults.csv", 'a')
	outfile.write("\n")
	outfile.close()

