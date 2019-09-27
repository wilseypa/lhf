/*
 * streamVR hpp + cpp extend the basePipe class for 
 * 
 */

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include <functional>
#include "slidingWindow.hpp"
#include "utils.hpp"
#include "readInput.hpp"


// basePipe constructor
slidingWindow::slidingWindow(){
	pipeType = "SlidingWindow";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket slidingWindow::runPipe(pipePacket inData){	
	utils ut;	
	readInput rp;
	
	// For this pipe, we construct a sub-pipeline:
	//		1. Read data vector by vector, push into slidingWindow evaluation
	//		2. IF the point is to be inserted, push into the simplexTree (LRU ordered)
	//		3. IF 100 points have passed, generate persistence intervals
	//
	// Loop this subpipeline until there's no more data
	
	std::vector<double> currentVector;
	rp.streamInit(inputFile);
	int pointCounter = 1;
	
	while(rp.streamRead(currentVector)){
		
		//Evaluate insertion into sliding window
		
		
		
		
		
		
		
		
		
		
		//Check if we've gone through 100 points
		if(pointCounter % 100 == 0){
			//Build and trigger remaining pipeline
			inData.complex->expandDimensions(dim);
			
			
			
		}
		pointCounter++;
		
	}
	//Probably want to trigger the remaining pipeline one last time...
	if((pointCounter - 1) % 100 != 0){
		
		inData.complex->expandDimensions(dim);
	}
	
	
	ut.writeLog("slidingWindow", "\tSuccessfully evaluated " + std::to_string(pointCounter) + " points");

	return inData;
}



// configPipe -> configure the function settings of this pipeline segment
bool slidingWindow::configPipe(std::map<std::string, std::string> configMap){
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = configMap["outputFile"].c_str();
	
	ut = utils(strDebug, outputFile);
	
	pipe = configMap.find("inputFile");
	if(pipe != configMap.end())
		inputFile = configMap["inputFile"].c_str();
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		epsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	else return false;
	
	ut.writeDebug("slidingWindow","Configured with parameters { input: " + configMap["inputFile"] + ", dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void slidingWindow::outputData(pipePacket inData){
	std::ofstream file ("output/" + pipeType + "_output.csv");
	
	if(inData.complex->simplexType == "simplexArrayList"){
		for(auto a : inData.complex->weightedGraph[1]){
			for(auto d : a.first){
				file << d << ",";
			}
			file << "\n";
		}
	}else{
		auto edges = inData.complex->getAllEdges(5);
		for (auto edge : edges){
			for (auto a : edge){
				for(auto d : a.first){
					file << d << ",";
				}
				file << a.second << "\n";
			}
		}
	}
	file << std::endl;
	
	file.close();
	return;
}
