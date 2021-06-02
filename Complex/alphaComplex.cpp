#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "alphaComplex.hpp"

template<typename nodeType>
alphaComplex<nodeType>::alphaComplex(double maxE, double maxD) : simplexArrayList<nodeType>::simplexArrayList(0,0) {
	this->simplexType = "alphaComplex";
	this->maxEpsilon = maxE;
	this->maxDimension = maxD;
}

template<typename nodeType>
alphaComplex<nodeType>::~alphaComplex(){
	this->simplexList.clear();
}

//Explicit Template Class Instantiation
template class alphaComplex<simplexNode>;
template class alphaComplex<alphaNode>;
template class alphaComplex<witnessNode>;

