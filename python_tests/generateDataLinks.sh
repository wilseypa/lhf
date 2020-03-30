#!/bin/bash

## README: This script is for HPC lab testing purposes only
##
##	Copy or move this script to the directory for tests to run
##		Run the script to link the currently tracked data sets
##		Run the multibatch script to start experiments on all data sets

## Use 
##		ln -s FILE LINK
## for each dataset to create symbolic links

## TODO: Generate this file?

if [ -f /testData/tda/static/UCI/camel.csv ]; then
	ln -s /testData/tda/static/UCI/camel.csv ./camel.csv
fi

if [ -f /testData/tda/static/UCI/lion.csv ]; then
	ln -s /testData/tda/static/UCI/lion.csv ./lion.csv
fi

if [ -f /testData/tda/static/UCI/flam.csv ]; then
	ln -s /testData/tda/static/UCI/flam.csv ./flam.csv
fi

if [ -f /testData/tda/static/UCI/elephant.csv ]; then
	ln -s /testData/tda/static/UCI/elephant.csv ./elephant.csv
fi

if [ -f /testData/tda/static/UCI/gesture_a1_va3.csv ]; then
	ln -s /testData/tda/static/UCI/gesture_a1_va3.csv ./gesture_a1_va3.csv
fi

if [ -f /testData/tda/static/Circles.csv ]; then
	ln -s /testData/tda/static/Circles.csv ./Circles.csv
fi

if [ -f /testData/tda/static/adHocData/twoMoons.csv ]; then
	ln -s /testData/tda/static/adHocData/twoMoons.csv ./twoMoons.csv
fi

if [ -f /testData/tda/static/adHocData/twoCircles.csv ]; then
	ln -s /testData/tda/static/adHocData/twoCircles.csv ./twoCircles.csv
fi

if [ -f /testData/tda/static/HAR/mergeda9leftlegMAG.csv ]; then
	ln -s /testData/tda/static/HAR/mergeda9leftlegMAG.csv ./mergeda9leftlegMAG.csv
fi

if [ -f /testData/tda/static/HAR/mergeda13leftlegMAG.csv ]; then
	ln -s /testData/tda/static/HAR/mergeda13leftlegMAG.csv ./mergeda13leftlegMAG.csv
fi

if [ -f /testData/tda/static/HAR/mergeda14leftlegMAG.csv ]; then
	ln -s /testData/tda/static/HAR/mergeda14leftlegMAG.csv ./mergeda14leftlegMAG.csv
fi

if [ -f /testData/tda/static/HAR/mergeda18leftlegMAG.csv ]; then
	ln -s /testData/tda/static/HAR/mergeda18leftlegMAG.csv ./mergeda18leftlegMAG.csv
fi

if [ -f /testData/tda/static/UCI/water-treatment.data_cleansed ]; then
	ln -s /testData/tda/static/UCI/water-treatment.data_cleansed ./water-treatment.data_cleansed
fi

if [ -f /testData/tda/static/UCI/iris-cleansed.data ]; then
	ln -s /testData/tda/static/UCI/iris-cleansed.data ./iris-cleansed.data
fi

if [ -f /testData/tda/static/otter_pointclouds/dragon_vrip.ply.txt_2000_.txt ]; then
	ln -s /testData/tda/static/otter_pointclouds/dragon_vrip.ply.txt_2000_.txt ./dragon_vrip.ply.txt_2000_.txt
fi

if [ -f /testData/tda/static/scikitTDAdata/torus.csv ]; then
	ln -s /testData/tda/static/scikitTDAdata/torus.csv ./torus.csv
fi

if [ -f /testData/tda/static/adHocData/4_rings.csv ]; then
	ln -s /testData/tda/static/adHocData/4_rings.csv ./4_rings.csv
fi

if [ -f /testData/tda/static/adHocData/spiral.csv ]; then
	ln -s /testData/tda/static/adHocData/spiral.csv ./spiral.csv
fi

if [ -f /testData/tda/static/scikitTDAdata/dsphere.csv ]; then
	ln -s /testData/tda/static/scikitTDAdata/dsphere.csv ./dsphere.csv
fi

if [ -f /testData/tda/static/scikitTDAdata/swiss_roll.csv ]; then
	ln -s /testData/tda/static/scikitTDAdata/swiss_roll.csv ./swiss_roll.csv
fi

if [ -f /testData/tda/static/otter_pointclouds/klein_bottle_pointcloud_new_900.txt ]; then
	ln -s /testData/tda/static/otter_pointclouds/klein_bottle_pointcloud_new_900.txt ./klein_bottle_pointcloud_new_900.txt
fi

if [ -f /testData/tda/static/scikitTDAdata/inf_sign.csv ]; then
	ln -s /testData/tda/static/scikitTDAdata/inf_sign.csv ./inf_sign.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/crosstrainMergedSegs45dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/crosstrainMergedSegs45dim.csv ./crosstrainMergedSegs45dim.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/jumpingMergedSegs45dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/jumpingMergedSegs45dim.csv ./jumpingMergedSegs45dim.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/stepperMergedSegs45dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/stepperMergedSegs45dim.csv ./stepperMergedSegs45dim.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/walkingMergedSegs45dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_45dim/walkingMergedSegs45dim.csv ./walkingMergedSegs45dim.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/crosstrainMergedSegs9dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/crosstrainMergedSegs9dim.csv ./crosstrainMergedSegs9dim.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/jumpingMergedSegs9dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/jumpingMergedSegs9dim.csv ./jumpingMergedSegs9dim.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/stepperMergedSegs9dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/stepperMergedSegs9dim.csv ./stepperMergedSegs9dim.csv
fi

if [ -f /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/walkingMergedSegs9dim.csv ]; then
	ln -s /testData/tda/static/chazalSubsample/subsamplePaperData/mergedData_9dim/walkingMergedSegs9dim.csv ./walkingMergedSegs9dim.csv
fi
