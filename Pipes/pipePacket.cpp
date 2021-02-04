#include <string>
#include <vector>
#include <typeinfo>
#include <iostream>
#include "pipePacket.hpp"

// pipePacket constructor, currently no needed information for the class constructor
pipePacket::pipePacket(const std::string& simplexType, const double epsilon, const int maxDim){
	std::map<std::string,std::string> blankConfig;
	blankConfig["dimensions"] = std::to_string(maxDim);
	blankConfig["epsilon"] = std::to_string(epsilon);
	
	if(complex != nullptr)
		delete complex;
	complex = simplexBase::newSimplex(simplexType, blankConfig);
}

pipePacket::pipePacket(std::map<std::string, std::string> configMap, const std::string& simplexType){
	
	if(complex != nullptr)
		delete complex;
	complex = simplexBase::newSimplex(simplexType, configMap);
}

std::string pipePacket::getStats(){
	std::string ret;
	ret += std::to_string(inputData.size()) + ",";
	ret += std::to_string(complex->simplexCount());
	
	return ret;
}

double pipePacket::getSize(){
	size_t size = 0;
	
	//1. Calculate size of original data
	for(auto row : workData){
		size += row.size() * sizeof(row[0]);
	}
	
	//2. Calculate size of input data
	for(auto row : inputData){
		size += row.size() * sizeof(row[0]);
	}
	
	//3. Calculate size of centroid labels
	size += centroidLabels.size() * sizeof(centroidLabels[0]);
	
	//4. Calculate size of the distance matrix
	for(auto row : distMatrix){
		size += row.size() * sizeof(row[0]);
	}
	
	//5. Calculate size of complex storage
	size += complex->getSize();
	
	//6. Calculate size of the boundaries
	for(auto row : boundaries){
		size += row.size() * sizeof(row.begin());
	}
	
	//7. Calculate size of weights
	size += weights.size() * sizeof(weights.begin());
	
	
	//8. Calculate size of bettiTable
	for(auto betti : bettiTable){
		size += betti.getSize();
	}
	
	//9. Calculate size of stats
	size += stats.size() * sizeof(std::string);
	
	//10. Calculate size of runLog
	size += runLog.size() * sizeof(std::string);
	
	return size;
}

