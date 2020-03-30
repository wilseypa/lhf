/*
 * 
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
#include "upscalePipe.hpp"
#include "utils.hpp"



// basePipe constructor
upscalePipe::upscalePipe(){
	pipeType = "Upscale";
	return;
}




// runPipe -> Run the configured functions of this pipeline segment
//
//	Upscale Pipe:
//			1. For each boundary vector identified
//				-Use index to find centroid label
//				-Upscale to original points from label
//				
//		
//
pipePacket upscalePipe::runPipe(pipePacket inData){
	
	std::vector<std::vector<std::vector<double>>> upscaledBoundaries;
	utils ut;
	
	//TODO: Update this to work with the bettiTable information
	//upscaledBoundaries = ut.separateBoundaryPartitions(wD->bettiTable);
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool upscalePipe::configPipe(std::map<std::string, std::string> configMap){
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
	
	configured = true;
	ut.writeDebug("upscale","Configured with parameters { debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

