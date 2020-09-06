/*
 * boundaryPipe hpp + cpp extend the basePipe class for calculating the 
 * minimum boundary from a distance matrix, used for upscaling the 
 * data for recomputation
 * 
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include "boundaryPipe.hpp"

// basePipe constructor
boundaryPipe::boundaryPipe(){
	pipeType = "Boundary";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//
pipePacket boundaryPipe::runPipe(pipePacket inData){
	
	// if(dim < 1)
	// 	return inData;
		
	// int count = 0;
		
	// std::vector<std::set<unsigned>> tempBetti;
	
	// std::cout << "Using scalar value: " << scalarV << std::endl;
	
	// //Extract the relevant boundaries from the centroid-approximated point cloud
	// for(auto bet : inData.bettiTable){
		
	// 	//If the betti is d1 or higher, or the distance is greater than the scalar
	// 	if(bet.bettiDim > 0 || bet.death - bet.birth >= scalarV){
			
			
	// 		std::set<unsigned> totalBoundary;
	// 		for(auto index : bet.boundary){
	// 			for(unsigned entry : index->simplex )
	// 				totalBoundary.insert(entry);
	// 		}
	// 		count++;
	// 		tempBetti.push_back(totalBoundary);
			
	// 	}
	// }
	
	// std::cout << "Boundaries updated: " << count << std::endl;
	
	// inData.boundaries = tempBetti;
	
	return inData;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void boundaryPipe::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_bettis_output.csv");
	
	file << inData.bettiOutput;
	
	file.close();
	return;
}

// configPipe -> configure the function settings of this pipeline segment
bool boundaryPipe::configPipe(std::map<std::string, std::string> configMap){
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
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end())
		dim = std::atoi(configMap["dimensions"].c_str());
	else return false;
	
	pipe = configMap.find("scalarV");
	if(pipe != configMap.end())
		scalarV = std::atof(configMap["scalarV"].c_str());
	else return false;
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	configured = true;
	ut.writeDebug("boundary","Configured with parameters { dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	
	return true;
}

