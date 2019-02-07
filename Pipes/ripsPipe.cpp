/*
 * ripsPipe hpp + cpp extend the basePipe class for calculating the 
 * VR complex from a neighborhood graph
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include "ripsPipe.hpp"
#include "utils.hpp"

// basePipe constructor
ripsPipe::ripsPipe(){
	pipeType = "Vietoris Rips (inductive)";
	return;
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket ripsPipe::runPipe(pipePacket inData){

	std::cout << "\tTest Complex Size: " << inData.workData.complex->edges.size() << std::endl;
	utils ut;
	for(auto z : inData.workData.complex->edges){
		ut.print1DVector(z);
	}
	
	std::cout << "EXPANDING DIMENSIONS" << std::endl;
	inData.workData.complex->expandDimensions(dim);
		
	std::cout << "\tComplex Size: " << inData.workData.complex->simplexCount() << std::endl;
	std::cout << "\tComplex Mem: " << inData.workData.complex->getSize() << std::endl;
	
	auto a = inData.workData.complex->getEdges(0,0);
	//std::cout << std::endl << std::endl;
	//for(auto z : a){
	//	std::cout << z.first << "\t";
	//	ut.print1DVector(z.second);
	//}
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool ripsPipe::configPipe(std::map<std::string, std::string> configMap){
	
	
	auto pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	
	return true;
}

