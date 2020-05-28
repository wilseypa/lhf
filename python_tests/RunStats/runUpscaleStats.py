import persim
import numpy as np
import sys
import time
import os

if not os.path.exists(outDir + '/upscaling'):
	os.makedirs(outDir + '/upscaling')
	
parser = argparse.ArgumentParser(description='Run iterative testing for TDA tools')
parser.add_argument('--SourcePers','-s',type=str, help='Baseline persistence intervals', default='Circles.csv')
parser.add_argument('--outDir', '-o', type=str, help='outDirectory', default='.')

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
print("Computing bottlenecks...")
bn = persim.bottleneck(originalPers, comparePers, matching=False)
bn_u = persim.bottleneck(originalPers, upscalePers, matching=False)

print("Computing heat kernel distance...")
h = persim.heat(originalPers, comparePers)
h_u = persim.heat(originalPers, upscalePers)

print("Computing wasserstein...")
ws = persim.wasserstein(originalPers, comparePers)
ws_u = persim.wasserstein(originalPers, upscalePers)
end = time.time()
#gh = persim.gromov_hausdorff(originalPers, comparePers)
#gh_u = persim.gromov_hausdorff(originalPers, upscalePers)

print(bn, bn_u, h, h_u, ws, ws_u)

stat_time = (end - start)

with open("upscaleStats.csv", 'a') as f:
	f.write(SourcePers + "," + outDir + ",".join(str(x) for x in [bn, bn_u, h, h_u, ws, ws_u, stat_time]) + "\n")
