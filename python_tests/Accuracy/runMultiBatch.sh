#!/bin/bash

#for file in camel.csv lion.csv flam.csv elephant.csv gesture_a1_va3.csv Circles.csv twoMoons.csv twoCircles.csv mergeda9leftlegMAG.csv mergeda13leftlegMAG.csv mergeda14leftlegMAG.csv mergeda18leftlegMAG.csv seeds_dataset_cleansed.txt water-treatment.data_cleansed iris-cleansed.data dragon_vrip.ply.txt_2000_.txt torus.csv 4_rings.csv spiral.csv dsphere.csv swiss_roll.csv klein_bottle_pointcloud_new_900.txt inf_sign.csv crosstrainMergedSegs45dim.csv jumpingMergedSegs45dim.csv stepperMergedSegs45dim.csv walkingMergedSegs45dim.csv crosstrainMergedSegs9dim.csv jumpingMergedSegs9dim.csv stepperMergedSegs9dim.csv walkingMergedSegs9dim.csv
for file in Circles.csv
#for file in 50_Points.csv
do
	for dim in 1 2 3 4 5
	do
		#Generate ground truths
		python runGroundTruths.py -f $file -e 2.0 -d $dim

		for cents in 25 50 75 100 125 150 175 200 250 500 1000
		do
			python runAccuracy.py -f $file -e 2.0 -d $dim -k $cents --np 5 -g ${file}_d${dim}_PH
		done
	done
done
