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
template<typename nodeType>
neighGraphPipe<nodeType>::neighGraphPipe(){
	this->pipeType = "neighGraph";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template<typename nodeType>
void neighGraphPipe<nodeType>::runPipe(pipePacket<nodeType> &inData){	
	
	//Iterate through each vector, inserting into simplex storage
	for(unsigned i = 0; i < inData.workData.size(); i++){
		if(!inData.workData[i].empty()){
			//insert data into the complex (SimplexArrayList, SimplexTree)
			inData.complex->insert();	
		}
	}

	return;
}

// configPipe -> configure the function settings of this pipeline segment
template<typename nodeType>
bool neighGraphPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();
	
	this->ut = utils(strDebug, this->outputFile);
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		this->epsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	else return false;
	
	this->configured = true;
	this->ut.writeDebug("neighGraphPipe","Configured with parameters { dim: " + std::to_string(dim) + " , eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
template<typename nodeType>
void neighGraphPipe<nodeType>::outputData(pipePacket<nodeType> &inData){
	std::ofstream file ("output/" + this->pipeType + "_output.csv");
	
	auto edges = inData.complex->getAllEdges();
	for (auto edge : edges){
		for (auto a : edge){
			for(auto d : a->simplex){
				file << d << ",";
			}
			file << a->weight << "\n";
		}
	}
	
	file << std::endl;
	
	file.close();
	return;
}

//Explicit Template Class Instantiation
template class neighGraphPipe<simplexNode>;
template class neighGraphPipe<alphaNode>;
template class neighGraphPipe<witnessNode>;
