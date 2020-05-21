import sys
import numpy as np
import persim
import argparse
import csv
import subprocess
import time
from sklearn import cluster as cs
import os
import string
np.set_printoptions(threshold=sys.maxsize) # for debugging
tdaLibraries = ['LHF','ripser']


''' Helper function to extract barcodes from Ripser output '''
def extractRipserData(barcodes, epsilon):
	retData = ""
	dim = -1

	for line in barcodes.split('\n'):

		if(line[0:4] == "pers"):
			dim += 1
		elif(line[0:2] == " [" and line[-2] == " "):
			retData += str(dim) + "," + line[2:-2] + str(epsilon) + "\n"
		elif(line[0:2] == " ["):
			retData += str(dim) + "," + line[2:-1] + "\n"

	return retData


''' MAIN LOOP HERE '''

''' Get CMD Arguments: [filename] [epsilon] [maxDim] [partitioner] [upscale 0/1] [centroids #] [reps]'''
parser = argparse.ArgumentParser(description='Run iterative testing for TDA tools')
parser.add_argument('--epsilon','-e',type=float, help='Max epsilon to compute PH up to', default=5)
parser.add_argument('--dim', '-d', type=int, help='Maximum homology dimension to compute (Hx)', default=1)
parser.add_argument('--filename', '-f', type=str, help='file to use', required=True)
parser.add_argument('--centroids', '-k', type=int, help='Number of centroids to compute', default=50)
parser.add_argument('--np','-np', type=int, help='Number of processors to use', default=1)
parser.add_argument('--reps','-r',type=int,help='Number of experiment repititions', default = 1)
parser.add_argument('--mode','-m',type=str,help='Mode to run LHF in', default='mpi')
parser.add_argument('--ground', '-g', type=str, help='Ground truths to compare against', required=True)
args = parser.parse_args()

epsilon = args.epsilon
maxDim = args.dim
centroids = args.centroids
reps = args.reps
mode = args.mode
noproc = args.np
filename = args.filename
comparePers = np.genfromtxt(args.ground, delimiter=',')
timings = {}

print "\nRunning " + mode + " experiments with:\nEpsilon: ", epsilon, "\nMaxDim", maxDim, "\nCentroids: ", centroids, "\nNProcs", noproc

''' Read input data '''
originalData = np.genfromtxt(filename, delimiter=',')
reducedData = originalData


'''Reduce the data if larger than max size'''
if(centroids < len(originalData)):
	
	reducedData = cs.KMeans(n_clusters=centroids, init='k-means++', n_jobs=-1).fit(originalData)
	originalLabels = reducedData.labels_
	reducedData = reducedData.cluster_centers_   

for i in range(0, int(reps)):

	''' Create folder for output '''
	outDir = filename + "_v"+str(centroids)+"_e"+str(epsilon) + "_d"+str(maxDim)+"/" + str(i) + "/"
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	if not os.path.isfile(os.getcwd() + "/aggResults.csv"):
		outfile = file(os.getcwd() + "/aggResults.csv", 'a')
		outfile.write('Filename,OutputPath,Vectors,H_Dim,PH_Library,Epsilon,PH_Time(s),BettiCount,BN,HK,WS,Stat_Time(s)\n')
		outfile.close()      
		
		
	''' Output all data to files '''
	np.savetxt(outDir + "originalData.csv", originalData, delimiter=',')
	np.savetxt(outDir + "processedData.csv", reducedData, delimiter=',')

	''' Run each TDA tool on the partitioned data '''
	for tdaType in tdaLibraries:
		
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
		outfile.write(filename + ',' + outDir[:-1] + "," + str(len(reducedData)) + "," + str(len(reducedData[0])) + "," + tdaType + "," + str(epsilon) + ",")
		outfile.close()

		''' RIPSER configuration and execution '''
		if(tdaType == "ripser"):
			if os.path.isfile("./ripser"):
					start = time.time()
					proc = subprocess.check_output(["./ripser", "--format", "point-cloud", "--dim", str(maxDim), "--threshold", str(epsilon), os.getcwd(
					)+"/"+outDir + "processedData.csv"])  # ,">>" + os.getcwd()+"/"+outDir+str(centroids)+"_Ripser.csv"])#, stdout=PIPE, stderr=PIPE)
					end = time.time()
					pers_time = (end - start)
					outdata = extractRipserData(proc, epsilon)
					outfile = file(os.getcwd() + "/" + outDir +
								   tdaType+"_Output.csv", 'w')
					outfile.write(outdata)
					outfile.close()
					
					outdata = np.genfromtxt(os.getcwd() + "/" + outDir +
								   tdaType+"_Output.csv", delimiter=',')

					start = time.time()
					print "Computing Metrics"
					bn = persim.bottleneck(outdata, comparePers, matching=False)

					h = persim.heat(outdata, comparePers)

					ws = persim.wasserstein(outdata, comparePers)
					end = time.time()
					stat_time=(end-start)
			

					outfile = file(os.getcwd() + "/aggResults.csv", 'a')
					outfile.write(str(pers_time) + "," + str(len(outdata)) + "," + str(bn) +"," + str(h) + "," + str(ws) + "," + str(stat_time) + ",")
					outfile.close()
			else:
				print "ripser processing failed"

		''' LHF configuration and execution '''
		if(tdaType == "LHF"):
			#try:
			start = time.time()
			print "mpiexec","-np",str(noproc),"./LHF", "-m",mode,"-d", str(maxDim), "-e", str(epsilon), "-k",str(centroids), "-s", "1", "-i", os.getcwd()+"/"+outDir + "originalData.csv"
			
			if(centroids == len(originalData)):
				proc = subprocess.check_output(["./LHF", "-m","fast","-d", str(maxDim+1), "-e", str(epsilon), "-s", "1", "-o", "output", "-i", os.getcwd()+"/"+outDir + "originalData.csv"])  
			else:
				proc = subprocess.check_output(["mpiexec","-np",str(noproc),"./LHF", "-m",mode,"-d", str(maxDim+1), "-e", str(epsilon), "-k",str(centroids), "-s", "1", "-o", "output", "-i", os.getcwd()+"/"+outDir + "originalData.csv"])  
			
			end = time.time()
			pers_time = (end - start)
			outdata = np.genfromtxt("./output.csv", delimiter=',')
			
			np.savetxt(os.getcwd() + "/" + outDir +tdaType+"_Output.csv", outdata, delimiter=',')
					
			
		   #except:
				#print "LHF processing failed"    
			
			start = time.time()
			print "Computing Metrics"
			bn = persim.bottleneck(outdata, comparePers, matching=False)

			h = persim.heat(outdata, comparePers)

			ws = persim.wasserstein(outdata, comparePers)
			end = time.time()
			stat_time=(end-start)
			
			outfile = file(os.getcwd() + "/aggResults.csv", 'a')
			outfile.write(str(pers_time) + "," + str(len(outdata)) + "," + str(bn) +"," + str(h) + "," + str(ws) + "," + str(stat_time) + ",")
			outfile.close()
			

		outfile = file(os.getcwd() + "/aggResults.csv", 'a')
		outfile.write("\n")
		outfile.close()
