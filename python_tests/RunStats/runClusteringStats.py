from sklearn import cluster as cs
from sklearn import metrics as met
import persim
import numpy as np
import sys
import time
import subprocess
import os
import argparse


#Upscaling - call Ripser (temporarily) using source data 
SourcePers = "Circles.csv"
outDir = "."
#need to change upscaling language to deal with different clustering algorithms
#need to find/add library to compute ARI, silhouette score, and Jaccard coeff
#to compare different clusterings
if not os.path.exists(outDir + '/clusterStats'):
	os.makedirs(outDir + '/clusterStats')
	

''' Get CMD Arguments: [SourcePers] [outDir] [epsilon] [dim] '''
parser = argparse.ArgumentParser(description='Run iterative testing for TDA tools')
parser.add_argument('--SourcePers','-s',type=str, help='Baseline persistence intervals', default='Circles.csv')
parser.add_argument('--outDir', '-o', type=str, help='outDirectory', default='.')
#.add_argument('--epsilon','-e',type=float, help='Max epsilon to compute PH up to', default=5)
#parser.add_argument('--dim', '-d', type=int, help='Maximum homology dimension to compute (Hx)', default=1)

args = parser.parse_args()

sourcePers = args.SourcePers
outDir = args.OutDir



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

baselineCluster = np.genfromtxt(outDir + "kmeans++/reducedData.csv")
print "Calculating ARI with k-means as the baseline comparison..."
#measure (dis)similarity between clusterings
#sklearn.metrics.adjusted_rand_score(labels_true, labels_pred)
#placeholder for now, need to confrim with nick
kmeansARI = met.adjusted_rand_score(baselineCluster, baselineCluster)
agglomWardARI = met.adjusted_rand_score(baselineCluster, outDir + "agglomerativeWard/reducedData.csv")
agglomSingleARi = met.adjusted_rand_score(baselineCluster, outDir + "agglomerativeSingle/reducedData.csv")
hdbscanARI = met.adjusted_rand_score(baselineCluster, outDir + "hdbscan/reducedData.csv")
randomARI = met.adjusted_rand_score(baselineCluster, outDir + "random/reducedData.csv")
print "Calculating Silhouette score " 
#silhouette sore - shows how close each point in one cluster is to points in neighboring clusters - used to find optimum num clusters
#sklearn.metrics.silhouette_score(X, labels, metric='euclidean', sample_size=None, random_state=None, **kwds)
#kmeansARI = met.silhouette_score(
#agglomWardARI = met.silhouette_score(
#agglomSingleARi = met.silhouette_score(
#hdbscanARI = met.silhouette_score(
#randomARI = met.silhouette_score(


print kmeansARI, agglomWardARI, agglomSingleARI, hdbscanARI, randomARI

print bn, bn_u, h, h_u, ws, ws_u

stat_time = (end - start)

with open("ClusteringStats.csv", 'a') as f:
	f.write(SourcePers + "," + outDir + ",".join(str(x) for x in [bn, bn_u, h, h_u, ws, ws_u, stat_time]) + "\n")
	f.write(kmeansARI + " kmeansARI\n", agglomWardARI + " agglomWardARI\n", agglomSingleARI + " agglomSingleARI\n", hdbscanARI + " hdbscanARI\n", randomARI +  " randomARI\n")
