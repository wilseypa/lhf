/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute  
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include "basePipe.hpp"
#include "distMatrixPipe.hpp"

// basePipe constructor
basePipe::basePipe(){
	

	return;
}

basePipe* basePipe::newPipe(const std::string &pipeT){
	std::cout << "Building pipeline: " << pipeT << std::endl;
	pipeType = pipeT;
	if(pipeType == "distMatrix"){
		std::cout << "test" << std::endl;
		return new distMatrixPipe();
	}
	
	return 0;
}

// runPipe -> Run the configured functions of this pipeline segment
std::vector<std::vector<double>> basePipe::runPipe(std::vector<std::vector<double>> inData){
	
	std::cout << "No run function defined for: " << pipeType << std::endl;
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool basePipe::configPipe(std::map<std::string, std::string> configMap){
	
	std::cout << "No configure function defined for: " << pipeType << std::endl;
	
	return true;
}

