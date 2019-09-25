import os
import sys
import csv
import subprocess
import numpy as np
import pandas as pd
np.set_printoptions(threshold=sys.maxsize)
#choose parameters to compare by using R, taken from cmd
#example comparisons:
#partitioner: RIPSER only: use kmeans++ as baseline: iterate through everything else..?
#partitioning size: RIPSER only: use max vector as baseline At each epsilon (param1) compare to all smaller vectors at corresponding Eps (param2)--- iterate by epsilon
#epsilon: RIPSER only: use max epsilon as baseline at each vector (param1) compare to all smaller epsilons at each vector (param2) ----iterate by vector
#library comparison: use Ripser at current vector and epsilon as baseline (param1) to compare to other libraries (param2) iterate by epsilon, then iterate by vector
#repition comparision: for each vector size, compare each repetition of each library --iterate by rep "/x"-- then iterate library, then iterate vector


#defaults
# Get CMD Arguments: [dataset] [library] [compareBy] [outputSeparate]
parentFile = "aggResults.csv"
dataset= "twoCircles"   
library= "ripser"   #Ripser, Eirene, GUDHI
dimension = 3
compareBy = "partitionSize"      # partitionSize, partitioner (preproc),  Epsilon, library, repetition, 
baseline = "Vectors" #baseline measurement. Set to Largest attribute of given parameter
outputSeparate = "Y"   #output separate file or append to aggResults.csv

if len(sys.argv) >=2:
	parentFile = sys.argv[1]
if len(sys.argv) >= 3:
	dataset = sys.argv[2]
if len(sys.argv) >= 4:
	library = sys.argv[3]
if len(sys.argv) >=5:
	dimension = sys.argv[4]
if len(sys.argv) >= 6:
	compareBy = sys.argv[5]
if len(sys.argv) >= 7:
	baseline = sys.argv[6]
if len(sys.argv) >= 8:
	outputSeparate = sys.argv[7]

path = os.getcwd()


#example path: aggResults.csv[OutputFolder]/0/EireneBettis        note: camel.csv_output_v10_e0.1/0 -> camel csv, vectors: 10, epsilon 0.1, repetition 0
if(library == "ripser"):
	resultsFile = "ripser_Output.csv"
elif(library == "Eirene"):
	resultsFile = "EireneBettis"

elif(library == "GUDHI"):
	resultsFile = "GUDHI_output.csv"

#path: fileName/0/  (EireneBettis)  (ripser_Output.csv) (GUDHI_Output.csv)
# Eirene -> cols 1,2 Ripser&Gudhi -> col 2,3



print "Running metrics with:\nDataSet: ", dataset, ":\nfrom parent file:", parentFile, "\ncomparing By:", compareBy, "\nLibrary:", library, "\nBaseline:", baseline, ":\noutput Separate?", outputSeparate

 
df = pd.read_csv(parentFile)     
#print df

df = df.dropna(axis='columns') #chopping off empty columns

#use panda dataframe to determine files to iterate through in the R script



#read cmd params to look for what to iterate through in aggResults to access corresponding file with Output data 


#Rscript scriptName Args

if(compareBy != "library" and baseline == "Vectors" ): #find effect of partitioning on PH
	 
	df = df[df.PH_Library == library]
	df = df[df.OutputPath.str.startswith(str(dataset))] # now only looking at chosen dataset
	#print df

	baselineVal  = max(df[str(baseline)]) #set baseline to largest value of current selected baseline parameter
	
	outfile = file(os.getcwd() + "/" + "Metrics_Output.csv", 'w')
	outfile.write("BaselineFile,CompareFile,Bottleneck0,Wasserstein0,Bottleneck1,Wasserstein1,Bottleneck2,Wasserstein2,Bottleneck3,Wasserstein3 \n")
	outfile.close()
	for i in range(0, len(df)):
		curEpsilon = df.iloc[i].Epsilon
		print curEpsilon
		
		try:
			baselineDiagram = df[(df[baseline] == baselineVal) & (df.Epsilon == curEpsilon) & (df.Preprocessor == "kmeans++")] #bandaid for now but fixes if theres no ripser output for 2000 vectors
		except:
			try:
				baselineVal = 1500
				baselineDiagram = df[(df[baseline] == baselineVal) & (df.Epsilon == curEpsilon) & (df.Preprocessor == "kmeans++")]
			except:
				try:
					baselineVal = 1000
					baselineDiagram = df[(df[baseline] == baselineVal) & (df.Epsilon == curEpsilon) & (df.Preprocessor == "kmeans++")]
				except:
					try:
						baselineVal = 750
						baselineDiagram = df[(df[baseline] == baselineVal) & (df.Epsilon == curEpsilon) & (df.Preprocessor == "kmeans++")]
					except:
						try:
							baselineVal = 500
							baselineDiagram = df[(df[baseline] == baselineVal) & (df.Epsilon == curEpsilon) & (df.Preprocessor == "kmeans++")]
						except:
							try:
								baselineVal = 250
								baselineDiagram = df[(df[baseline] == baselineVal) & (df.Epsilon == curEpsilon) & (df.Preprocessor == "kmeans++")]
							except:
								print "no baseline vector found"

		#print baselineDiagram.iloc[0].Preprocessor
		baselineFolder = baselineDiagram.iloc[0].OutputPath #+ baselineDiagram.iloc[0].Preprocessor   #["OutputPath"]
		compareFolder = df.iloc[i].OutputPath
		#print baselineFolder
		print compareFolder
		baselineFile = baselineFolder + "/" + resultsFile
		print baselineFile
		compareFile = compareFolder + "/" + resultsFile
		outfile = file(os.getcwd() + "/" + "Metrics_Output.csv", 'a')
		outdata = subprocess.check_output(["Rscript", "TDA_metrics.r", str(baselineFile), str(compareFile),  dimension, str(compareBy)]) #store output from R script as string
		writeData = (baselineFile + "," + compareFile + "," + outdata + "\n")
		outfile.write(writeData)
		outfile.close()

	#	print (library + "library R script failed")
elif(compareBy != "library" and baseline == "Epsilon" ):  #find effect of epsilon on PH
	try: 
		df = df[df.TDALibrary == library]
		df = df[df.OutputPath.str.startswith(dataset)] # now only looking at chosen dataset
		baselineVal  = max(df[baseline]) #set baseline to largest value of current selected baseline parameter
		
		for i in range(0, len(df)):
			curVector = df.iloc[i].Vectors
			baselineDiagram = df[(df[baseline] == baselineVal) & (df.Vectors == curVector)] 
			baselineFolder = baselineDiagram["OutputPath"]
			subprocess.check_output(["Rscript", "TDA_metrics.r",  baselineFolder, compareFolder, resultsFile, dimension, compareBy]) #store output from R script as string
	except:
		print (library + "library R script failed")
else:  #if comparing libraries.... baseline should be ripser
	try:
		df = df[df.OutputFolder.str.startswith(dataset)]
		subprocess.check_output(["Rscript", "TDA_metrics.r", parentFile, compareBy]) 
	except:
		print "Library compare script failed"
	

