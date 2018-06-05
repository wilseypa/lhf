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
#include "neighGraphPipe.hpp"
#include "bettiPipe.hpp"
#include "ripsPipe.hpp"

// basePipe constructor
basePipe::basePipe(){
	

	return;
}

basePipe* basePipe::newPipe(const std::string &pipeT){
	std::cout << "Building pipeline: " << pipeT << std::endl << std::endl;
	pipeType = pipeT;
	if(pipeType == "distMatrix"){
		return new distMatrixPipe();
	} else if (pipeType == "neighGraph"){
		return new neighGraphPipe();
	} else if (pipeType == "betti"){
		return new bettiPipe();
	} else if (pipeType == "rips"){
		return new ripsPipe();
	}
	
	return 0;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket basePipe::runPipe(pipePacket inData){
	
	std::cout << "No run function defined for: " << pipeType << std::endl;
	
	return inData;
}	

// configPipe -> configure the function settings of this pipeline segment
bool basePipe::configPipe(std::map<std::string, std::string> configMap){
	
	std::cout << "No configure function defined for: " << pipeType << std::endl;
	
	return true;
}

