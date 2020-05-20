#!/bin/bash

#for file in camel.csv lion.csv flam.csv elephant.csv gesture_a1_va3.csv Circles.csv twoMoons.csv twoCircles.csv mergeda9leftlegMAG.csv mergeda13leftlegMAG.csv mergeda14leftlegMAG.csv mergeda18leftlegMAG.csv seeds_dataset_cleansed.txt water-treatment.data_cleansed iris-cleansed.data dragon_vrip.ply.txt_2000_.txt torus.csv 4_rings.csv spiral.csv dsphere.csv swiss_roll.csv klein_bottle_pointcloud_new_900.txt inf_sign.csv crosstrainMergedSegs45dim.csv jumpingMergedSegs45dim.csv stepperMergedSegs45dim.csv walkingMergedSegs45dim.csv crosstrainMergedSegs9dim.csv jumpingMergedSegs9dim.csv stepperMergedSegs9dim.csv walkingMergedSegs9dim.csv

for dim in 1 2 3 4 5
do
	for procs in 1 2 3 4 #5 6 7 8
	do
		for centroids in 10 20 30 40 50 60 70 80 90 100
		do
			for points in 100 250 500 750 1000
			do
				python runScalability.py -e 2.0 -d $dim -n $points -k $centroids --np $procs
			done
		done
	done
done
