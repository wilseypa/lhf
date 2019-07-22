#include <string>
#include <vector>
#include <typeinfo>
#include <iostream>
#include "pipePacket.hpp"

// pipePacket constructor, currently no needed information for the class constructor
pipePacket::pipePacket(const std::string& simplexType, const double epsilon, const int maxDim){
	simplexBase *bs = new simplexBase(epsilon, maxDim);
	complex = bs->newSimplex(simplexType);
}

double pipePacket::getSize(){
	
	//Calculate size of original data
	
	size_t size = 0;
	
	for(auto row : originalData){
		for(auto index : row)
			size += sizeof(index);
	}	
	
	//Calculate size of complex storage
	size += complex->getSize();
	
	return size;
}
