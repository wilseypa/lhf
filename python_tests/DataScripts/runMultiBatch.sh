#!/bin/bash

export PYTHONPATH='/home/ubuntu/git/GUDHI/cython/'

for file in Circles.csv
do
	#Start by generating standardized/normalized data sets

	python scaleColumns.py $file

	#Next generate the histogram and data plots

	python dataPlots.py --filename $file
	python dataPlots.py --filename $file-standardized
	python dataPlots.py --filename $file-scaled

	#Finally attempt PH on the full dataset

	python tdaPlots.py --filename $file
	python tdaPlots.py --filename $file-standardized
	python tdaPlots.py --filename $file-scaled

done


