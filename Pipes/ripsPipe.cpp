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
template<typename nodeType>
ripsPipe<nodeType>::ripsPipe(){
	this->pipeType = "ripsPipe";
	return;
}


// runPipe -> Run the configured functions of this pipeline segment
template<typename nodeType>
void ripsPipe<nodeType>::runPipe(pipePacket<nodeType> &inData){
	
	inData.complex->expandDimensions(dim);
		
	this->ut.writeDebug("ripsPipe","Expanded Complex Size: " + std::to_string(inData.complex->simplexCount()));
	this->ut.writeDebug("ripsPipe", "Expanded Complex Mem: " + std::to_string(inData.complex->getSize()));
	
	/*if(collapse == "true" || collapse == "1"){
		inData.complex->reduceComplex();
		
		this->ut.writeDebug("ripsPipe","Reduced Complex Size: " + std::to_string(inData.complex->simplexCount()));
		this->ut.writeDebug("ripsPipe", "Reduced Complex Mem: " + std::to_string(inData.complex->getSize()));
	}*/
	return;
}


// configPipe -> configure the function settings of this pipeline segment
template<typename nodeType>
bool ripsPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
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
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	pipe = configMap.find("collapse");
	if(pipe != configMap.end())
		this->collapse = configMap["collapse"];
	
	this->configured = true;
	this->ut.writeDebug("ripsPipe","Configured with parameters { dim: " + std::to_string(dim) + " , debug: " + strDebug + ", outputFile: " + this->outputFile + ", collapse: " + this->collapse + " }");
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
template<typename nodeType>
void ripsPipe<nodeType>::outputData(pipePacket<nodeType> &inData){
	std::ofstream file;
	
	if(inData.complex->simplexType == "simplexArrayList"){
		file.open("output/" + this->pipeType + "_output.csv");
		for (int i = 0; i < inData.complex->simplexList.size(); i++){
			for(auto a : inData.complex->simplexList[i]){
				for(auto d : a->simplex){
					file << d << ",";
				}
				file << a->weight << "\n";
			}
		}
		
		file.close();
	}
	return;
}

template class ripsPipe<simplexNode>;
template class ripsPipe<alphaNode>;
