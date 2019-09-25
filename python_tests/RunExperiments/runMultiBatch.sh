#!/bin/bash

export PYTHONPATH='/home/ubuntu/git/GUDHI/cython/'

for file in camel.csv lion.csv flam.csv elephant.csv gesture_a1_va3.csv twoCircles.csv twoMoons.csv Circles.csv
do
	for points in 10 25 50 75 100 250 500 750 1000 1500 2000
	do
		for edge in 0.1 0.25 0.5 1.0
		do
			python runExperiments.py $file $edge 10 kmeans++ 0 $points 5
		done
	done
done


