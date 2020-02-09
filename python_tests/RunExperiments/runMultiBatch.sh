#!/bin/bash

export PYTHONPATH='/home/ubuntu/git/GUDHI/cython/'

for file in Circles.csv
do
	for points in 10 25 50 75 100 #250 500 750 1000
	do
		for edge in 1.0
		do
			python runExperiments.py -f $file -e $edge -d 2 -k $points
		done
	done
done


