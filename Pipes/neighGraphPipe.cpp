/*
 * neighGraphPipe hpp + cpp extend the basePipe class for building the 
 * neighborhood graph from data input in the complex structure configured
 * 
 */

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include <functional>
#include "neighGraphPipe.hpp"


// basePipe constructor
neighGraphPipe::neighGraphPipe(){
	pipeType = "neighGraph";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket neighGraphPipe::runPipe(pipePacket inData){	
	
	//Iterate through each vector, inserting into simplex storage
	for(unsigned i = 0; i < inData.originalData.size(); i++){
		if(!inData.originalData[i].empty()){
			//insert data into the complex (SimplexArrayList, SimplexTree)
			inData.complex->insert(inData.originalData[i]);	
		}
	}

	return inData;
}

// configPipe -> configure the function settings of this pipeline segment
bool neighGraphPipe::configPipe(std::map<std::string, std::string> configMap){
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
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		epsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	else return false;
	
	configured = true;
	ut.writeDebug("neighGraphPipe","Configured with parameters { dim: " + std::to_string(dim) + " , eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void neighGraphPipe::outputData(pipePacket inData){
	std::ofstream file ("output/" + pipeType + "_output.csv");
	
	if(inData.complex->simplexType == "simplexArrayList"){
		if(inData.complex->weightedGraph.size() > 1){
			for(auto a : inData.complex->weightedGraph[1]){
				for(auto d : a.first){
					file << d << ",";
				}
				file << a.second << "\n";
			}
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
