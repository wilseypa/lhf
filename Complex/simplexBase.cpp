#include <string>
#include <vector>
#include <iostream>
#include <typeinfo>
#include "simplexBase.hpp"
#include "simplexTree.hpp"
#include "simplexArrayList.hpp"
#include "indSimplexTree.hpp"

simplexBase::simplexBase(){return;}
simplexBase::simplexBase(double maxE, int maxDim){
	maxEpsilon = maxE;	
	maxDimension = maxDim;	
}

void simplexBase::setDistanceMatrix(std::vector<std::vector<double>> _distMatrix){
	distMatrix = _distMatrix;
	return;
}

// simplexTree constructor, currently no needed information for the class constructor
simplexBase* simplexBase::newSimplex(const std::string &simplexT){
	simplexType = simplexT;
	
	if(simplexType == "simplexTree"){
		return new simplexTree(maxEpsilon, distMatrix, maxDimension);
	} else if (simplexType == "simplexArrayList"){	
		return new simplexArrayList(maxEpsilon, distMatrix);
	} else if (simplexType == "indSimplexTree"){
		return new indSimplexTree(maxEpsilon, distMatrix, maxDimension);
	}
	return 0;
}


std::vector<std::vector<unsigned>> simplexBase::getEdges(int dim, double epsilon){
	std::cout << "No get edges function defined for: " << simplexType << std::endl;
	std::vector<std::vector<unsigned>> a;
	return a;
}

std::vector<std::vector<std::pair<std::set<unsigned>,double>>> simplexBase::getAllEdges(double epsilon){
	std::cout << "No get edges function defined for: " << simplexType << std::endl;
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> a;
	return a;
}

double simplexBase::getSize(){
	std::cout << "No size function defined for: " << simplexType << std::endl;
	return 0;
}

void simplexBase::insert(std::vector<double>&){
	std::cout << "No insert function defined for: " << simplexType << std::endl;
	return;
}

void simplexBase::find(std::vector<double>){
	std::cout << "No find function defined for: " << simplexType << std::endl;
	return;
}

int simplexBase::vertexCount(){
	std::cout << "No vertexCount function defined for: " << simplexType << std::endl;
	return -1;
}

int simplexBase::simplexCount(){
	std::cout << "No simplexCount function defined for: " << simplexType << std::endl;
	return -1;
}

void simplexBase::outputSimplex(){
	std::cout << "No outputSimplex function defined for: " << simplexType << std::endl;
	return;
}

void simplexBase::expandDimensions(int dim){
	std::cout << "No expandDimensions function defined for: " << simplexType << std::endl;
	return;
}
	
