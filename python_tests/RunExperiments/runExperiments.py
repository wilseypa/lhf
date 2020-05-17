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

def geometricMean(tempData, originalData):
    
    tracker = max(tempData.labels_)
    sumSize = tracker +1
    dim = len(originalData[0])
    summedClusters = np.zeros([ sumSize, dim])
    mean = np.zeros([sumSize])
    for t in range(tracker+1):
        for i in range(len(originalData)):
            #print t
            if(tempData.labels_[i] == t):
                #print summedClusters
                summedClusters[t] += originalData[i]


    unique, counts = np.unique(tempData.labels_, return_counts=True)
    countsDict = dict(zip(unique, counts))
    countsTemp = np.array(counts)
    countsTemp2 = np.array([ [counts]*dim ]).T
    countsTrick = np.reshape(countsTemp2, (sumSize, dim)) #juggling the stupid numpy arrays so we can trick the vectors for geo means
    
    
    mean = summedClusters/countsTrick

    reducedData = mean
    
    return reducedData
def geometricMeanDenoise(tempData, originalData): #identical to geometric mean but for algorithms that have noise 
    
    tracker = max(tempData.labels_)
    sumSize = tracker +1
    dim = len(originalData[0])
    summedClusters = np.zeros([ sumSize, dim])
    mean = np.zeros([sumSize])
    for t in range(tracker+1):
        for i in range(len(originalData)):
            #print t
            if(tempData.labels_[i] == t):
                #print summedClusters
                summedClusters[t] += originalData[i]


    unique, counts = np.unique(tempData.labels_, return_counts=True)
    countsDict = dict(zip(unique, counts))
    countsTemp = np.array(counts[1:]) #start at 1 because -1 is noise points
    countsTemp2 = np.array([ [countsTemp]*dim ]).T
    countsTrick = np.reshape(countsTemp2, (sumSize, dim)) 
    
    
    mean = summedClusters/countsTrick

    reducedData = mean
    
    return reducedData    
    
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
    return reducedData


''' Helper function to extract barcodes, boundaries from Eirene output '''

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


''' Helper function to extract barcodes from Ripser output '''

def extractRipserData(barcodes, epsilon, outfile):
    epsilon = str(epsilon)
    bettiCount = 0
    for dim in range(len(barcodes)):
        bettiCount += len(barcodes[dim])
        for line in barcodes[dim]:
            outfile.write(str(dim) + "," + str(line[0]) + "," + (epsilon if line[1] == np.inf else str(line[1])) + "\n")

    return bettiCount


''' Helper function for writing GUDHI persistence to file '''


def outputGudhiPersistence(data, label, maxedge):
    buf = ""

    for i in data:
        tmp = str(i).replace('(', '')
        tmp = str(tmp).replace(' ', '')
        tmp = str(tmp).replace('inf', maxedge)
        buf += tmp[:-2] + "\n"

    out = open(outDir+'GUDHI_Output.csv', 'w')
    out.write(buf)
    out.close()
    return


''' MAIN LOOP HERE '''

''' Get CMD Arguments: [filename] [epsilon] [maxDim] [partitioner] [upscale 0/1] [centroids #] [reps]'''
parser = argparse.ArgumentParser(description='Run iterative testing for TDA tools')
parser.add_argument('--filename','-f',type=str, help='Point cloud filename', default='Circles.csv')
parser.add_argument('--partitioner', '-p', type=str, help='Partitioner for point cloud reduction', default='kmeans++')
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
timings = {}
originalLabels = []
neighSize = args.neighSize

print("Running experiments with:\nFilename: ", fileName, "\nEpsilon: ", epsilon, "\nMaxDim", maxDim, "\nPartitioner: ",  partitioner, "\nUpscale: ", upscale, "\nCentroids: ", centroids, "\nMinCluster HDBSCAN", neighSize)


''' Read input data '''
originalData = np.genfromtxt(fileName, delimiter=',')
reducedData = originalData


for i in range(0, int(reps)):

    ''' Create folder for output '''
    outDir = os.path.basename(fileName) + "_output_v"+str(centroids)+"_e"+str(epsilon) + "/" + str(i) + "/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        
    if not os.path.isfile(os.getcwd() + "/aggResults.csv"):
        outfile = open(os.getcwd() + "/aggResults.csv", 'a')
        outfile.write('Filename,OutputPath,Vectors,Dimensions,Preprocessor,PreprocessingTime(s),PH_Library,Epsilon,PH_Time(s),BettiCount\n')
        outfile.close()      

    timings[partitioner]=0;
    ''' Partition data and store centroids, labels '''
    
    if(len(originalData[0]) < maxDim):
        maxDim = len(originalData[0])
    

    if(centroids < len(originalData)):
        start = time.time()
    
        if(partitioner == "kmeans++"): #baseline partitioner
            reducedData = cs.KMeans(n_clusters=centroids, init='k-means++', n_jobs=-1).fit(originalData)
            originalLabels = reducedData.labels_
            reducedData = reducedData.cluster_centers_          
        elif(partitioner == "agglomerativeWard"):  #ward linkage --> minimizes the variance of the clusters being merged
            tempData = cs.AgglomerativeClustering(n_clusters=centroids, linkage="ward").fit(originalData)
            reducedData = geometricMean(tempData, originalData)
            #print reducedData
        elif(partitioner == "agglomerativeSingle"):  #single linkage --> the minimum of the distances between all observations of the two sets
            tempData = cs.AgglomerativeClustering(n_clusters=centroids, linkage="single").fit(originalData)
            reducedData = geometricMean(tempData, originalData)
            #print reducedData
            #HDBSCAN - updated version of DBSCAN that iterates through all epsilons and builds a hierarchy of density... tune  minClusterSize to 
            #change size of data reduction ie 3 min Pts = low Reduction 100 min Pts = high Reduction
        elif(partitioner == "hdbscan"): 
            clusterer = hdbscan.HDBSCAN(algorithm='best', alpha=1.0, approx_min_span_tree=True, gen_min_span_tree=False, leaf_size=40, metric='euclidean', min_cluster_size=neighSize, min_samples=None, p=None)
            tempData = clusterer.fit(originalData)  
            reducedData = geometricMeanDenoise(tempData, originalData)
        elif(partitioner == "random"):
            reducedData = np.array(random.sample(originalData, centroids)) #random.sample = no duplicates
            #print len(reducedData)
            
        end = time.time()

        timings[partitioner] = (end - start)


    ''' Output all data to files '''
    np.savetxt(outDir + "reducedData.csv", reducedData, delimiter=',')
    np.savetxt(outDir + "reducedDataPivot.csv", reducedData.transpose(), delimiter=',')
    np.savetxt(outDir + "originalLabels.csv", originalLabels, delimiter=',')


    ''' Run each TDA tool on the partitioned data '''
    for tdaType in tdaLibraries:
        
        #Check if the last line was returned, keeping our data clean
        isRet = True;
        with open(os.getcwd() + "/aggResults.csv", 'rb+') as f:
            f.seek(-1,2)
            a = f.read()
            if(a is not "\n"):
                isRet = False;
        
        outfile = open(os.getcwd() + "/aggResults.csv", 'a')
        if(not isRet):
            outfile.write('')
        outfile.write(fileName + ',' + outDir[:-1] + "," + str(len(reducedData)) + "," + str(len(reducedData[0])) + "," + str(timings[partitioner]) + "," + str(partitioner) + "," + tdaType + "," + str(epsilon) + ",")
        outfile.close()

        ''' Eirene configuration and execution '''
        if(tdaType == "Eirene"):
            try:
                proc = subprocess.check_output(
				["julia", "runEirene.jl", os.getcwd() + "/" + outDir, str(maxDim), str(epsilon)])
                pers_time = 0
                outBarcodes, outBoundaries = extractEireneData(proc, epsilon)            

                outfile = open(os.getcwd() + "/" + outDir + tdaType + "_Boundaries.txt","w")
                outBoundaries = str(outBoundaries).replace('),','\n')
                outfile.write(outBoundaries.translate(None,'[]()set'))
                outfile.close()            
                	
                outfile = open(os.getcwd() + "/" + outDir + tdaType + "_Output.csv", "w")
                outfile.write(outBarcodes)
                outfile.close()

                with open(outDir + '/EireneTime') as f:
                	pers_time = f.read()
                	
                outfile = open(os.getcwd() + "/aggResults.csv", 'a')
                outfile.write(str(pers_time) + "," + str(str(outBarcodes.count("\n"))) + ",")
                outfile.close()
            except:
                print("Eirene processing failed")

        ''' RIPSER configuration and execution '''
        if(tdaType == "ripser"):
            try:
                start = time.time()
                barcodes = ripser(reducedData, maxdim = maxDim-1, thresh = epsilon)['dgms']
                end = time.time()
                pers_time = (end - start)

                outfile = open(os.getcwd() + "/" + outDir + tdaType+"_Output.csv", 'w')
                bettiCount = extractRipserData(barcodes, epsilon, outfile)
                outfile.close()

                outfile = open(os.getcwd() + "/aggResults.csv", 'a')
                outfile.write(str(pers_time) + "," + str(bettiCount) + ",")
                outfile.close()
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
                pers_time = (end - start)
                outputGudhiPersistence(pers, str(centroids) + "_", str(epsilon))

                outfile = open(os.getcwd() + "/aggResults.csv", 'a')
                outfile.write(str(pers_time) + "," + str(len(pers)) + ",")
                outfile.close()
            except:
                print("GUDHI processing failed")
                
                
        ''' LHF configuration and execution '''
        if(tdaType == "LHF"):
            #try:
            start = time.time()
            proc = subprocess.check_output(["../../LHFmain/LHF", "-m","fast","-d", str(maxDim+1), "-e", str(epsilon), "-i", os.getcwd()+"/"+outDir + "reducedData.csv"])  


            end = time.time()
            pers_time = (end - start)
            outdata = np.genfromtxt("./output.csv.csv")
            
            outfile = open(os.getcwd() + "/" + outDir +tdaType+"_Output.csv", 'w')
                    
            np.savetxt(os.getcwd() + "/" + outDir +tdaType+"_Output.csv", outdata)
            # outfile.write(outdata)
                    
            outfile.close()
                    
            outfile = open(os.getcwd() + "/aggResults.csv", 'a')
            outfile.write(str(pers_time) + "," + str(outdata.size - 1) + ",")
            outfile.close()
           #except:
                #print "LHF processing failed"    


        outfile = open(os.getcwd() + "/aggResults.csv", 'a')
        outfile.write("\n")
        outfile.close()

