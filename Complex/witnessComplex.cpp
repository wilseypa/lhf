#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "witnessComplex.hpp"

template<typename nodeType>
witnessComplex<nodeType>::witnessComplex(double maxE, double maxD) : simplexArrayList<nodeType>::simplexArrayList(0,0) {
	this->simplexType = "alphaComplex";
	this->maxEpsilon = maxE;
	this->maxDimension = maxD;
}

template<typename nodeType>
witnessComplex<nodeType>::~witnessComplex(){
	this->simplexList.clear();
}
