/*
 * ripsPipe hpp + cpp extend the basePipe class for calculating the 
 * VR complex from a neighborhood graph
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
#include "ripsPipe.hpp"
#include "utils.hpp"

// basePipe constructor
ripsPipe::ripsPipe(){
	pipeType = "Rips_Inductive";
	return;
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket ripsPipe::runPipe(pipePacket inData){
	utils ut;
	
	inData.workData.complex->expandDimensions(dim);
		
	ut.writeLog("ripsPipe","Expanded Complex Size: " + std::to_string(inData.workData.complex->simplexCount()));
	ut.writeLog("ripsPipe", "Expanded Complex Mem: " + std::to_string(inData.workData.complex->getSize()));
	
	inData.workData.complex->reduceComplex();
	
	ut.writeLog("ripsPipe","Reduced Complex Size: " + std::to_string(inData.workData.complex->simplexCount()));
	ut.writeLog("ripsPipe", "Reduced Complex Mem: " + std::to_string(inData.workData.complex->getSize()));
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool ripsPipe::configPipe(std::map<std::string, std::string> configMap){
	
	auto pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	ut.writeDebug("ripsPipe","Configured with parameters { dim: " + std::to_string(dim) + "}");
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void ripsPipe::outputData(pipePacket inData){
	std::ofstream file;
	
	if(inData.workData.complex->simplexType == "simplexArrayList"){
		file.open("output/" + pipeType + "_output.csv");
		for (int i = 0; i < inData.workData.complex->weightedGraph.size(); i++){
			for(auto a : inData.workData.complex->weightedGraph[i]){
				for(auto d : a){
					file << d << ",";
				}
				file << "\n";
			}
		}
		
		file.close();
	}
	return;
}

