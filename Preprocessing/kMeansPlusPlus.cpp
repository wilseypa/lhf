/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute  
 * 
 */

#include <chrono>
#include <string>
#include <iostream>
#include <vector>
#include "kMeansPlusPlus.hpp"

// basePipe constructor
kMeansPlusPlus::kMeansPlusPlus(){
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket kMeansPlusPlus::runPreprocessor(pipePacket inData){
	
	std::cout << "No run function defined for: " << procName << std::endl;
	
	return inData;
}	

// configPipe -> configure the function settings of this pipeline segment
bool kMeansPlusPlus::configPreprocessor(std::map<std::string, std::string> configMap){
	
	std::cout << "No configure function defined for: " << procName << std::endl;
	
	return true;
}

