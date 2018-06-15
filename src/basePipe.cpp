/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute  
 * 
 */

#include <chrono>
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

// runPipeWrapper -> wrapper for timing of runPipe and other misc. functions
pipePacket basePipe::runPipeWrapper(pipePacket inData){
	pipePacket retData;
	
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	retData = runPipe(inData);
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();
	
	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
	
	//Output the time and memory used for this pipeline segment
	std::cout << "Pipeline " << pipeType << " executed in " << (elapsed.count()/1000.0) << " seconds (physical time)" << std::endl;
	
	auto dataSize = retData.getSize();
	auto unit = "Bytes";
	
	//Convert large datatypes (GB, MB, KB)
	if(dataSize > 1000000000){
		//Convert to GB
		dataSize = dataSize/1000000000;
		unit = "GigaBytes";
	} else if(dataSize > 1000000){
		//Convert to MB
		dataSize = dataSize/1000000;
		unit = "MegaBytes";
	} else if (dataSize > 1000){
		//Convert to KB
		dataSize = dataSize/1000;
		unit = "KiloBytes";
	}
	
	std::cout << "\tData size: " << dataSize << " " << unit << std::endl << std::endl;
	
	
	return retData;
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

