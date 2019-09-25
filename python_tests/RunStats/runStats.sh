#!/bin/bash

for dir in `find . -mindepth 2 -maxdepth 2 -type d | grep camel-standardized`
do
	echo $dir
	echo "${dir%%_*}"
	python runUpscaleStats.py "camel-standardized.csv_output_v25_e1.0/0/ripser_Output.csv" $dir 3.0 3
done
 
