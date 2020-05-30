from sklearn import cluster as cs
import sys
import numpy as np
import argparse
import csv
import subprocess
import time
import os
import gudhi
import random
import hdbscan
import string
from ripser import ripser
np.set_printoptions(threshold=sys.maxsize) # for debugging
tdaLibraries = ['LHF','ripser','GUDHI']



## HELPER FUNCTIONS --------------------------------------------#

'''Geometric mean of clusters'''
def clusterCenters(data, labels):
    n = max(labels)
    sums = np.zeros([n+1, len(data[0])])
    counts = np.zeros([n+1])
    for i in range(len(data)):
        if(labels[i] != -1):
            sums[labels[i]] += data[i]
            counts[labels[i]] += 1

    return sums/counts[:, None]

'''
def corePointExtract(tempData, originalData, centroids):
    unique, counts = np.unique(tempData.labels_, return_counts=True)
    countsDict = dict(zip(unique, counts))
    #print countsDict
    clusterIndex = max(tempData.labels_)
    clusters = [ [] for cluster in range(clusterIndex+1)]
   
    
    for c in range(0, clusterIndex+1):
        for i in range(len(originalData)):
            if (tempData.labels_[i] == c):
                clusters[c].append(originalData[i])
    rawClusters = np.asarray(clusters)

    centroidFrac = centroids/(clusterIndex + 1)
    #print centroidFrac
    tempClusters = []
    for c in range(0, clusterIndex +1):
        tempClusters.append(cs.KMeans(n_clusters = centroidFrac, init="k-means++", n_jobs = -1).fit(rawClusters[c]).cluster_centers_)

    finalClusters = np.vstack(tempClusters)  
  
    reducedData = finalClusters
    return reducedData'''


''' Helper function to extract barcodes, boundaries from Eirene output '''
'''
def extractEireneData(procOutput, epsilon):
    retBarcodes = ""
    retBoundaries = []
    dim = -1
    
    #Start with parsing barcodes and cycles out
    tempProc = proc.split('Barcodes')[1]
    Barcodes = tempProc.split('Cycles')[0]
    Cycles = tempProc.split('Cycles')[1]
    
    #print "Cycles: ",Cycles
    
    #Parse and format the barcodes
    for line in Barcodes.split('\n'):
        dim = line.split(':')[0]
        if(line.find("Array") < 0 and len(line) > 10):
            rem = line.split('[')[1]
            rem = rem.replace("; ","\n")
            rem = rem.replace(" ",",")
            rem = rem.replace("Inf",str(epsilon))
            for subline in rem.split("\n"):
                retBarcodes += dim + "," + subline+ "\n"
        
        
    retBarcodes = retBarcodes.replace("]","");

    #Parse and format the boundaries
    for line in Cycles.split('\n'):
        if(len(line) > 5):
            dim = line.split(':')[0]
            for subline in line.split(':[')[1].split(';'):
                a = set()
                for index in subline.split(' '):
                    index = index.replace(']','')
                    if(len(index) > 0):
                        a.add(int(index))
                
                boundary = subline.split(", ");
                retBoundaries.append(a)
                
    curBounds = []
    #Reduce boundary sets
    for i in range(len(retBoundaries)):
        insert = True
        for j in range(len(curBounds)):
            if(len(retBoundaries[i].intersection(curBounds[j])) > 0):
                curBounds[j] = curBounds[j].union(retBoundaries[i])
                insert = False
                break
        if(insert):
            curBounds.append(retBoundaries[i])
    
    print("RetBounds:\n",retBoundaries)
    retBoundaries = curBounds
    
    print("CurBounds:\n",retBoundaries)
    return retBarcodes, retBoundaries
'''

''' Helper function to extract barcodes from Ripser output '''
def extractRipserData(barcodes, epsilon, outfile):
    bettiCount = 0
    for dim in range(len(barcodes)):
        bettiCount += len(barcodes[dim])
        for line in barcodes[dim]:
            outfile.write(str(dim) + "," + str(line[0]) + "," + (epsilon if line[1] == np.inf else str(line[1])) + "\n")

    return bettiCount


''' Helper function for writing GUDHI persistence to file '''
def outputGudhiPersistence(data, epsilon, outfile):
    for line in data:
        outfile.write(str(line[0]) + ',' + str(line[1][0]) + ',' + (epsilon if line[1][1] == np.inf else str(line[1][1])) + "\n")


''' MAIN LOOP HERE '''

''' Get CMD Arguments: [filename] [epsilon] [maxDim] [partitioner] [upscale 0/1] [centroids #] [reps]'''
parser = argparse.ArgumentParser(description='Run iterative testing for TDA tools')
parser.add_argument('--filename','-f',type=str, help='Point cloud filename', default='Circles.csv')
parser.add_argument('--partitioner', '-p', type=str, help='Partitioner for point cloud reduction', default=None)
parser.add_argument('--epsilon','-e',type=float, help='Max epsilon to compute PH up to', default=5)
parser.add_argument('--dim', '-d', type=int, help='Maximum homology dimension to compute (Hx)', default=1)
parser.add_argument('--upscale', '-u', type=bool, help='Set whether to upscale data during processing', default=0)
parser.add_argument('--centroids', '-k', type=int, help='Number of centroids to compute', default=50)
parser.add_argument('--reps','-r',type=int,help='Number of experiment repititions', default = 1)
parser.add_argument('--neighSize', '-n', type=int, help='Minimum cluster size for HDBSCAN', default = 3)
args = parser.parse_args()

fileName = args.filename
partitioner = args.partitioner
epsilon = args.epsilon
maxDim = args.dim
upscale = args.upscale
centroids = args.centroids
reps = args.reps
originalLabels = []
neighSize = args.neighSize

print("Running experiments with:\nFilename: ", fileName, "\nEpsilon: ", epsilon, "\nMaxDim", maxDim, "\nPartitioner: ",  partitioner, "\nUpscale: ", upscale, "\nCentroids: ", centroids, "\nMinCluster HDBSCAN", neighSize)


''' Read input data '''
originalData = np.genfromtxt(fileName, delimiter=',')
reducedData = originalData

if not os.path.isfile(os.getcwd() + "/aggResults.csv"):
    outfile = open(os.getcwd() + "/aggResults.csv", 'a')
    outfile.write('Filename,OutputPath,Vectors,Dimensions,Preprocessor,PreprocessingTime,PH_Library,Epsilon,PH_Time,BettiCount\n')
    outfile.close()

for i in range(reps):

    ''' Create folder for output '''
    outDir = os.path.basename(fileName) + "_output_v"+str(centroids)+"_e"+str(epsilon) + "/" + str(i) + "/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    partitionTime = 0
    ''' Partition data and store centroids, labels '''
    
    if(len(originalData[0]) < maxDim):
        maxDim = len(originalData[0])
    

    if(partitioner != None):
        start = time.time()
    
        if(partitioner == "kmeans++"): #baseline partitioner
            reducedData = cs.KMeans(n_clusters=centroids, init='k-means++', n_jobs=-1).fit(originalData)
            originalLabels = reducedData.labels_
            reducedData = reducedData.cluster_centers_          
        elif(partitioner == "agglomerativeWard"):  #ward linkage --> minimizes the variance of the clusters being merged
            labels = cs.AgglomerativeClustering(n_clusters=centroids, linkage="ward").fit_predict(originalData)
            reducedData = clusterCenters(originalData, labels)
        elif(partitioner == "agglomerativeSingle"):  #single linkage --> the minimum of the distances between all observations of the two sets
            labels = cs.AgglomerativeClustering(n_clusters=centroids, linkage="single").fit_predict(originalData)
            reducedData = clusterCenters(originalData, labels)
            #HDBSCAN - updated version of DBSCAN that iterates through all epsilons and builds a hierarchy of density... tune  minClusterSize to 
            #change size of data reduction ie 3 min Pts = low Reduction 100 min Pts = high Reduction
        elif(partitioner == "hdbscan"): 
            clusterer = hdbscan.HDBSCAN(algorithm='best', alpha=1.0, approx_min_span_tree=True, gen_min_span_tree=False, leaf_size=40, metric='euclidean', min_cluster_size=neighSize, min_samples=None, p=None)
            labels = clusterer.fit_predict(originalData)  
            reducedData = clusterCenters(originalData, labels)
        elif(partitioner == "random"):
            reducedData = originalData[np.random.choice(originalData.shape[0], centroids, replace=False), :]
            
        end = time.time()

        partitionTime = (end - start)


    ''' Output all data to files '''
    np.savetxt(outDir + "reducedData.csv", reducedData, delimiter=',')
    np.savetxt(outDir + "reducedDataPivot.csv", reducedData.transpose(), delimiter=',')
    np.savetxt(outDir + "originalLabels.csv", originalLabels, delimiter=',')


    ''' Run each TDA tool on the partitioned data '''
    for tdaType in tdaLibraries:
        persTime = 0
        bettiCount = 0
        file = os.getcwd() + "/" + outDir +tdaType+"_Output.csv"

        ''' Eirene configuration and execution '''
        '''
        if(tdaType == "Eirene"):
            try:
                proc = subprocess.check_output(
				["julia", "runEirene.jl", os.getcwd() + "/" + outDir, str(maxDim), str(epsilon)])
                persTime = 0
                outBarcodes, outBoundaries = extractEireneData(proc, epsilon)            

                outfile = open(os.getcwd() + "/" + outDir + tdaType + "_Boundaries.txt","w")
                outBoundaries = str(outBoundaries).replace('),','\n')
                outfile.write(outBoundaries.translate(None,'[]()set'))
                outfile.close()            
                	
                outfile = open(os.getcwd() + "/" + outDir + tdaType + "_Output.csv", "w")
                outfile.write(outBarcodes)
                outfile.close()

                with open(outDir + '/EireneTime') as f:
                	persTime = f.read()
                	
                bettiCount = outBarcodes.count("\n")
            except:
                print("Eirene processing failed") '''

        ''' RIPSER configuration and execution '''
        if(tdaType == "ripser"):
            try:
                start = time.time()
                barcodes = ripser(reducedData, maxdim = maxDim-1, thresh = epsilon)['dgms']
                end = time.time()
                persTime = (end - start)

                with open(file, 'w') as outfile:
                    bettiCount = extractRipserData(barcodes, str(epsilon), outfile)
            except:
                print("ripser processing failed")

        ''' GUDHI configuration and execution '''
        if(tdaType == "GUDHI"):
            try:
                start = time.time()
                rips = gudhi.RipsComplex(points=reducedData, max_edge_length=epsilon)
                simplex_tree = rips.create_simplex_tree(max_dimension=maxDim)
                pers = simplex_tree.persistence()
                end = time.time()
                persTime = (end - start)
                bettiCount = len(pers)

                with open(file, 'w') as outfile:
                    outputGudhiPersistence(pers, str(epsilon), outfile)
            except:
                print("GUDHI processing failed")
                
                
        ''' LHF configuration and execution '''
        if(tdaType == "LHF"):
            try:
                start = time.time()
                proc = subprocess.check_output(["../../LHFmain/LHF", "-m","fast","-d", str(maxDim), "-e", str(epsilon), "-i", os.getcwd()+"/"+outDir + "reducedData.csv"])  
                end = time.time()
                persTime = (end - start)

                outdata = np.genfromtxt("./output.csv.csv", delimiter=',')                
                np.savetxt(file, outdata, delimiter=',')
                bettiCount = len(outdata)
            except:
                print("LHF processing failed")  


        with open(os.getcwd() + "/aggResults.csv", 'a') as outfile:
            outfile.write(fileName + ',' + outDir[:-1] + "," + str(len(reducedData)) + "," + str(len(reducedData[0])) + "," + str(partitioner) + "," + str(partitionTime) + "," + tdaType + "," + str(epsilon) + "," + str(persTime) + "," + str(bettiCount) + ",\n")
