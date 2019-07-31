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
	pipeType = "ripsPipe";
	return;
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket ripsPipe::runPipe(pipePacket inData){
	
	inData.complex->expandDimensions(dim);
		
	ut.writeDebug("ripsPipe","Expanded Complex Size: " + std::to_string(inData.complex->simplexCount()));
	ut.writeDebug("ripsPipe", "Expanded Complex Mem: " + std::to_string(inData.complex->getSize()));
	
	if(collapse == "true" || collapse == "1"){
		inData.complex->reduceComplex();
		
		ut.writeDebug("ripsPipe","Reduced Complex Size: " + std::to_string(inData.complex->simplexCount()));
		ut.writeDebug("ripsPipe", "Reduced Complex Mem: " + std::to_string(inData.complex->getSize()));
	}
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool ripsPipe::configPipe(std::map<std::string, std::string> configMap){
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
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	pipe = configMap.find("collapse");
	if(pipe != configMap.end())
		collapse = configMap["collapse"];
	
	ut.writeDebug("ripsPipe","Configured with parameters { dim: " + std::to_string(dim) + " , debug: " + strDebug + ", outputFile: " + outputFile + ", collapse: " + collapse + " }");
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void ripsPipe::outputData(pipePacket inData){
	std::ofstream file;
	
	if(inData.complex->simplexType == "simplexArrayList"){
		file.open("output/" + pipeType + "_output.csv");
		for (int i = 0; i < inData.complex->weightedGraph.size(); i++){
			for(auto a : inData.complex->weightedGraph[i]){
				for(auto d : a.first){
					file << d << ",";
				}
				file << a.second << "\n";
			}
		}
		
		file.close();
	}
	return;
}

