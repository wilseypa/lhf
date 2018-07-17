#include <string>
#include <vector>
#include <iostream>
#include "simplexBase.hpp"
#include "simplexTree.hpp"
#include "simplexArrayList.hpp"

simplexBase::simplexBase(){
	
}

// simplexTree constructor, currently no needed information for the class constructor
simplexBase* simplexBase::newSimplex(const std::string &simplexT){
	simplexType = simplexT;
	
	std::cout << "Creating new simplex structure: " << simplexT << std::endl;
	
	if(simplexType == "simplexTree"){
		return new simplexTree();
	} else if (simplexType == "simplexArrayList"){
		return new simplexArrayList();
	}
	return 0;
}

double simplexBase::getSize(){
	std::cout << "No size function defined for: " << simplexType << std::endl;
	return 0;
}

void simplexBase::insert(std::vector<double>){
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
	
