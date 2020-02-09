
## Normalize the column values of a input .csv file using functions from sklearn.preprocessing.

## Two types of normalization are supported, namely:
##     MinMaxScaler: scale the data to the range [-1.0, 1.0]
##     StandardScaler: set the values so that the mean is 0.0 and the std deviation is 1.0.
##                     This is also called z-score normalization.
## The input values are assumed to be numeric but the sklearn functions are setup so that NaN will have no
## impact on the computation and will  be left in place.

## This script is very simple and has no real sanity checking.  It will read the filename off the command
## line, strip the last 4 characters from the filename (assuming them to be '.csv') and use this string to
## write two files: one with the suffix '-scaled.csv' containing the MinMaxScaler modified data, and one with
## the suffix '-standardized.csv' containing the StandardScaler modified data.

import sys
import numpy as np
import sklearn.preprocessing
import argparse
from argparse import RawTextHelpFormatter

####----------------------------------------------------------------------------------------------------
#### let us begin
####----------------------------------------------------------------------------------------------------

# process the arguments on the command line
argparser = argparse.ArgumentParser(description='''Normalize the features (columns) of the input file using sklearn.preprocessing.  Two types of normalization are performed, namely:

  MinMaxScaler: Scale the columns values of a data file to a range of [-1.0, 1.0].
  StandardScaler (known as z-score normalization): Set column values such that the mean = 0.0 and the std deviation = 1.0.

The output files will be written with suffixes: -scaled.csv, -standardized.csv''',
                                    formatter_class=RawTextHelpFormatter)
argparser.add_argument('fileName', help='Name of input file to scale.')

args = argparser.parse_args()

inFile = args.fileName
if inFile is None :
    print 'Missing input filename....aborting'
    sys.exit()

rawData = np.loadtxt(inFile, delimiter = ",", comments="#")

# scale each feature (column of data) to a range of (-1.0,1.0)
scaler = sklearn.preprocessing.MinMaxScaler(feature_range=(-1.0,1.0))
scaledData = scaler.fit_transform(rawData)

outFile = inFile[0:len(inFile)-4] + "-scaled.csv"

with open(outFile, 'w') as of :
    np.savetxt(of, scaledData, delimiter=',')
of.close()

standardizer = sklearn.preprocessing.StandardScaler().fit(rawData)
standardizedData = standardizer.transform(rawData)

outFile = inFile[0:len(inFile)-4] + "-standardized.csv"

with open(outFile, 'w') as of :
    np.savetxt(of, standardizedData, delimiter=',')
of.close()

