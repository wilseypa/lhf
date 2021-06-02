/*
 * distMatrix hpp + cpp extend the basePipe class for calculating the 
 * distance matrix from data input
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
#include <set>
#include <algorithm>
#include "distMatrixPipe.hpp"
#include "utils.hpp"

// basePipe constructor
template<typename nodeType>
distMatrixPipe<nodeType>::distMatrixPipe(){
	this->pipeType = "DistMatrix";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template<typename nodeType>
void distMatrixPipe<nodeType>::runPipe(pipePacket<nodeType> &inData){
	//Store our distance matrix
	if(inData.distMatrix.size() > 0) inData.distMatrix.clear();
	inData.distMatrix.resize(inData.workData.size(), std::vector<double>(inData.workData.size(),0));
	//Iterate through each vector, create lower
	for(unsigned i = 0; i < inData.workData.size(); i++){
		//Grab a second vector to compare to 
		for(unsigned j = i+1; j < inData.workData.size(); j++){
			//Calculate vector distance 
			inData.distMatrix[i][j] = this->ut.vectors_distance(inData.workData[i],inData.workData[j]);
		}
	}

	for(unsigned i = 0; i < inData.workData.size(); i++){
		double r_i = 0;
		for(unsigned j = 0; j < inData.workData.size(); j++) r_i = std::max(r_i, inData.distMatrix[std::min(i, j)][std::max(i, j)]);
		enclosingRadius = std::min(enclosingRadius, r_i);
	}

	inData.complex->setDistanceMatrix(&inData.distMatrix);
	inData.complex->setEnclosingRadius(enclosingRadius);

	this->ut.writeDebug("distMatrix", "\tDist Matrix Size: " + std::to_string(inData.distMatrix.size()) + " x " + std::to_string(inData.distMatrix.size()));
	return;
}


// configPipe -> configure the function settings of this pipeline segment
template<typename nodeType>
bool distMatrixPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
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
		this->enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else return false;

	this->configured = true;
	this->ut.writeDebug("distMatrixPipe","Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");
	
	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template<typename nodeType>
void distMatrixPipe<nodeType>::outputData(pipePacket<nodeType> &inData){
	std::ofstream file;
	file.open("output/" + this->pipeType + "_output.csv");
	
	for(std::vector<double> a : inData.distMatrix){
		for(auto d : a){
			file << d << ",";
		}
		file << "\n";
	}
	
	file.close();
	return;
}

//Explicit Template Class Instantiation
template class distMatrixPipe<simplexNode>;
template class distMatrixPipe<alphaNode>;
