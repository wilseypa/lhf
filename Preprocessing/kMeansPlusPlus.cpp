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
	//Arguments - How many clusters? Convergence value
	
	
	//Initialize centroids (Plus plus mechanism with kmeans - Hartigan, Wong)
	
	
	//Iterate over the data, compute WCSSE -> Minimize WCSSE
	
	
	
	//Convergence, return the centroids, the cluster assignments
	
	
	
	
	
	std::cout << "No run function defined for: " << procName << std::endl;
	
	return inData;
}	

// configPipe -> configure the function settings of this pipeline segment
bool kMeansPlusPlus::configPreprocessor(std::map<std::string, std::string> configMap){
	
	std::cout << "No configure function defined for: " << procName << std::endl;
	
	return true;
}

