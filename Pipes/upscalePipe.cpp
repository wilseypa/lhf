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
	
	std::vector<std::vector<double>> currentData;
		
	// Iterate each set of boundary points
	for(auto bound : inData.boundaries) {
		
		for(auto point : bound){
			
			currentData.push_back(inData.originalData[inData.originalLabels[point]]);
		}
	}
	inData.reducedData = currentData;
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool upscalePipe::configPipe(std::map<std::string, std::string> configMap){
	//auto pipe = configMap.find("dimensions");
	//if(pipe != configMap.end())
	//	dim = std::atoi(configMap["dimensions"].c_str());
	//else return false;
	
	return true;
}

