#!/usr/bin/python3

## a collection of curve manipulation utilities for our tda explorations.

import sys
import os
import numpy as np


def subtractCurve(orig, sub):
	"""
	Subtract two unaligned curves (ECC / Betti) where the first column represents the filtration
		value (weight) and the second column is the measured value
		
    Parameters
    ----------
	orig : numpy.array(n, 2)
		original curve, where [-1] represents weight and [0] represents the value
		
	sub : numpy.array(m, 2)
		curve to subtract from original, where [-1] represents weight and [0] represents
			the value		
		
	Return
    ------
    numpy.array(p, 2)
		Result of subtraction, (orig - sub). Number of results, p < 
	"""
	output = []
	subIndex = 0
	curSub = 0;
	
	for i in range(len(orig)):
		
		tout = orig[i].copy()
		
		#Iterate through orig; check if cur weight >= sub weight
		if tout[-1] >= sub[subIndex][-1]:
			tout[0] = tout[0] - sub[subIndex][0]
			curSub = sub[subIndex][0]
			subIndex += 1
		#interpolate until next change
		else:
			tout[0] -= sub[subIndex][0]
	
		output.append(tout)
		
	return output
	
def addCurve(orig, sub):
	"""
	Add two unaligned curves (ECC / Betti) where the first column represents the filtration
		value (weight) and the second column is the measured value
		
    Parameters
    ----------
	orig : numpy.array(n, 2)
		original curve, where [-1] represents weight and [0] represents the value
		
	sub : numpy.array(m, 2)
		curve to subtract from original, where [-1] represents weight and [0] represents
			the value		
		
	Return
    ------
    numpy.array(p, 2)
		Result of addition, (orig + sub)
	"""
	
	output = []
	subIndex = 0
	curSub = 0;
	
	for i in range(len(orig)):
		
		tout = orig[i].copy()
		
		#Iterate through orig; check if cur weight >= sub weight
		if tout[-1] >= sub[subIndex][-1]:
			tout[0] = tout[0] + sub[subIndex][0]
			curSub = sub[subIndex][0]
			subIndex += 1
		#interpolate until next change
		else:
			tout[0] += sub[subIndex][0]
	
		output.append(tout)
		
	return output	
