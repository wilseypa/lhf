import persim
import numpy as np
import sys
import time
import subprocess
import os


#Upscaling - call Ripser (temporarily) using source data 
SourcePers = "Circles.csv"
outDir = "."
#need to change upscaling language to deal with different clustering algorithms
#need to find/add library to compute ARI, silhouette score, and Jaccard coeff
#to compare different clusterings
if not os.path.exists(outDir + '/upscaling'):
	os.makedirs(outDir + '/upscaling')
	
if len(sys.argv) >= 2:
    SourcePers = sys.argv[1]
if len(sys.argv) >= 3:
	outDir = sys.argv[2]
if len(sys.argv) >= 4:
    epsilon = float(sys.argv[3])
if len(sys.argv) >= 5:
	dim = int(sys.argv[4])

originalPers = np.genfromtxt(SourcePers, delimiter=',')
comparePers = np.genfromtxt(outDir + "/Eirene_Output.csv", delimiter=',')
try:
	upscalePers = np.genfromtxt(outDir + "/upscaledPersistence.csv", delimiter = ',')
except:
	upscalePers = np.genfromtxt(outDir + "/ripser_Output.csv", delimiter = ',')
	

start = time.time()
print "Computing bottlenecks..."
bn = persim.bottleneck(originalPers, comparePers, matching=False)
bn_u = persim.bottleneck(originalPers, upscalePers, matching=False)

print "Computing heat kernel distance..."
h = persim.heat(originalPers, comparePers)
h_u = persim.heat(originalPers, upscalePers)

print "Computing wasserstein..."
ws = persim.wasserstein(originalPers, comparePers)
ws_u = persim.wasserstein(originalPers, upscalePers)
end = time.time()
#gh = persim.gromov_hausdorff(originalPers, comparePers)
#gh_u = persim.gromov_hausdorff(originalPers, upscalePers)

print bn, bn_u, h, h_u, ws, ws_u

stat_time = (end - start)

with open("upscaleStats.csv", 'a') as f:
	f.write(SourcePers + "," + outDir + ",".join(str(x) for x in [bn, bn_u, h, h_u, ws, ws_u, stat_time]) + "\n")
