/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute  
 * 
 */

#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "basePipe.hpp"
#include "distMatrixPipe.hpp"
#include "neighGraphPipe.hpp"
#include "bettiPipe.hpp"
#include "ripsPipe.hpp"
#include "upscalePipe.hpp"
#include "boundaryPipe.hpp"

// basePipe constructor
basePipe::basePipe(){
	return;
}

basePipe* basePipe::newPipe(const std::string &pipeT){
	std::cout << "Building pipeline: " << pipeT << std::endl;
	pipeType = pipeT;
	if(pipeType == "distMatrix"){
		return new distMatrixPipe();
	} else if (pipeType == "neighGraph"){
		return new neighGraphPipe();
	} else if (pipeType == "betti"){
		return new bettiPipe();
	} else if (pipeType == "rips"){
		return new ripsPipe();
	} else if (pipeType == "upscale"){
		return new upscalePipe();
	} else if (pipeType == "boundary"){
		return new boundaryPipe();
	}
	
	return 0;
}

// runPipeWrapper -> wrapper for timing of runPipe and other misc. functions
pipePacket basePipe::runPipeWrapper(pipePacket inData){
	
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	inData = runPipe(inData);
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();
	
	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
	
	//Output the time and memory used for this pipeline segment
	std::cout << "\tPipeline " << pipeType << " executed in " << (elapsed.count()/1000.0) << " seconds (physical time)" << std::endl;
	
	auto dataSize = inData.getSize();
	auto unit = "B";
	
	//Convert large datatypes (GB, MB, KB)
	if(dataSize > 1000000000){
		//Convert to GB
		dataSize = dataSize/1000000000;
		unit = "GB";
	} else if(dataSize > 1000000){
		//Convert to MB
		dataSize = dataSize/1000000;
		unit = "MB";
	} else if (dataSize > 1000){
		//Convert to KB
		dataSize = dataSize/1000;
		unit = "KB";
	}
	
	std::cout << "\t\tData size: " << dataSize << " " << unit << std::endl << std::endl;
	
	inData.stats += pipeType + "," + std::to_string(elapsed.count()/1000.0) + "," + std::to_string(dataSize) + "," + unit + "," + std::to_string(inData.workData.complex->vertexCount()) + "," + std::to_string(inData.workData.complex->simplexCount()) + "\n";
	
	
	std::cout << "ret..." << std::endl;
	if(debug){
		std::cout << "ret...2" << std::endl;
		outputData(inData);
	}
	
	return inData;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void basePipe::outputData(pipePacket inData){
	std::cout << "No output function defined for: " << pipeType << std::endl;
	
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");
	
	for (auto a : inData.workData.originalData){
		for (auto d : a){
			file << std::to_string(d) << ",";
		}
		file << "\n";
	}
	
	file.close();
	return;
}
	


// runPipe -> Run the configured functions of this pipeline segment
pipePacket basePipe::runPipe(pipePacket inData){
	
	std::cout << "No run function defined for: " << pipeType << std::endl;
	
	return inData;
}	

// configPipe -> configure the function settings of this pipeline segment
bool basePipe::configPipe(std::map<std::string, std::string> configMap){
	
	std::cout << "No configure function defined for: " << pipeType << std::endl;

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	else return false;
	
	return true;
}

