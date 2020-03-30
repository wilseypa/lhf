#!/bin/bash

export PYTHONPATH='/home/ubuntu/git/GUDHI/cython/'

for file in camel.csv lion.csv flam.csv elephant.csv gesture_a1_va3.csv Circles.csv twoMoons.csv twoCircles.csv mergeda9leftlegMAG.csv mergeda13leftlegMAG.csv mergeda14leftlegMAG.csv mergeda18leftlegMAG.csv seeds_dataset_cleansed.txt water-treatment.data_cleansed iris-cleansed.data dragon_vrip.ply.txt_2000_.txt torus.csv 4_rings.csv spiral.csv dsphere.csv swiss_roll.csv klein_bottle_pointcloud_new_900.txt inf_sign.csv crosstrainMergedSegs45dim.csv jumpingMergedSegs45dim.csv stepperMergedSegs45dim.csv walkingMergedSegs45dim.csv crosstrainMergedSegs9dim.csv jumpingMergedSegs9dim.csv stepperMergedSegs9dim.csv walkingMergedSegs9dim.csv
do
	for points in 10 25 50 75 100 250 500 750 1000
	do
		for edge in 1.0
		do
			python runExperiments.py -f $file -e $edge -d 2 -k $points
		done
	done
done


