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
#include "ripsPipe.hpp"
#include "upscalePipe.hpp"
#include "persistencePairs.hpp"
#include "slidingWindow.hpp"
#include "fastPersistence.hpp"
#include "incrementalPersistence.hpp"
#include "naiveWindow.hpp"

basePipe* basePipe::newPipe(const std::string &pipeType, const std::string &complexType){
	utils ut;
	ut.writeDebug("basePipe","Building pipeline: " + pipeType + " for " + complexType);

	if(pipeType == "distMatrix"){
		return new distMatrixPipe();
	} else if (pipeType == "neighGraph"){
		return new neighGraphPipe();
	} else if (pipeType == "rips"){
		return new ripsPipe();
	} else if (pipeType == "upscale"){
		std::cout << "Building upscale" << std::endl;
		return new upscalePipe();
	} else if (pipeType == "persistence"){
		return new persistencePairs();
	} else if (pipeType == "slidingwindow" || pipeType == "sliding"){
		return new slidingWindow();
	} else if (pipeType == "fastPersistence" || pipeType == "fast"){
		return new fastPersistence();
	} else if (pipeType == "incrementalPersistence" || pipeType == "inc"){
		return new incrementalPersistence();
	} else if (pipeType == "naivewindow" || pipeType == "naive"){
		return new naiveWindow();
	}
	
	return 0;
}

// runPipeWrapper -> wrapper for timing of runPipe and other misc. functions
void basePipe::runPipeWrapper(pipePacket &inData){
	
	//Check if the pipe has been configured
	if(!configured){
		ut.writeLog(pipeType,"Pipe not configured");
		std::cout << "Pipe not configured" << std::endl;
		return;
	}
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	runPipe(inData);
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();
	
	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
	
	//Output the time and transient memory used for this pipeline segment
	ut.writeLog(pipeType,"\tPipeline " + pipeType + " executed in " + std::to_string(elapsed.count()/1000.0) + " seconds (physical time)");
	
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
	
	inData.stats += pipeType + "," + std::to_string(elapsed.count()/1000.0) + "," + std::to_string(dataSize) + "," + unit + "," + std::to_string(inData.complex->vertexCount()) + "," + std::to_string(inData.complex->simplexCount()) + "\n";
	
	std::string ds = std::to_string(dataSize);
	ut.writeLog(pipeType,"\t\tData size: " + ds + " " + unit + "\n");
	
	if(debug)
		outputData(inData);
	
	return;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void basePipe::outputData(pipePacket &inData){
	ut.writeDebug("basePipe","No output function defined for: " + pipeType);
	
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");
	
	for (auto a : inData.workData){
		for (auto d : a){
			file << std::to_string(d) << ",";
		}
		file << "\n";
	}
	
	file.close();
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
void basePipe::runPipe(pipePacket &inData){
	ut.writeError("basePipe","No run function defined for: " + pipeType);
	
	return;
}	

// configPipe -> configure the function settings of this pipeline segment
bool basePipe::configPipe(std::map<std::string, std::string> &configMap){
	ut.writeDebug("basePipe","No configure function defined for: " + pipeType);

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	
	return true;
}
