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

template<typename nodeType>
distMatrixPipe<nodeType>::distMatrixPipe(){
    /**
	    distMatrixPipe()
	 
		@brief Class constructor
		@tparam nodeType The data type of the simplex node.
	*/
	this->pipeType = "DistMatrix";
	return;
}

template<typename nodeType>
void distMatrixPipe<nodeType>::runPipe(pipePacket<nodeType> &inData){
	/**
	    runPipe(pipePacket<nodeType> &inData)
	 
		@brief Constructs the distance matrix from the point cloud. Uses parameters set in respective configPipe.
		@tparam nodeType The data type of the simplex node.
		@param inData The pipePacket data being used in the pipeline.
	*/
	if(inData.distMatrix.size() > 0) 
		inData.distMatrix.clear();
	inData.distMatrix.resize(inData.workData.size(), std::vector<double>(inData.workData.size(),0));
	
	for(unsigned i = 0; i < inData.workData.size(); i++){
		for(unsigned j = i+1; j < inData.workData.size(); j++){
			//Vector distance between i, j
			inData.distMatrix[i][j] = this->ut.vectors_distance(inData.workData[i],inData.workData[j]);
		}
	}

	for(unsigned i = 0; i < inData.workData.size(); i++){
		double r_i = 0;
		
		for(unsigned j = 0; j < inData.workData.size(); j++) 
			r_i = std::max(r_i, inData.distMatrix[std::min(i, j)][std::max(i, j)]);
			
		enclosingRadius = std::min(enclosingRadius, r_i);
	}
	
	if(inData.complex->complexType == "betaComplex" && (this->betaMode == "lune" || this->betaMode == "circle"))
		inData.incidenceMatrix = this->ut.betaNeighbors(inData.inputData,beta,betaMode);

	inData.complex->setDistanceMatrix(&inData.distMatrix);
	inData.complex->setEnclosingRadius(enclosingRadius);
	inData.complex->setIncidenceMatrix(&inData.incidenceMatrix);

	this->ut.writeDebug("distMatrix", "\tDist Matrix Size: " + std::to_string(inData.distMatrix.size()) + " x " + std::to_string(inData.distMatrix.size()));
	return;
}


// configPipe -> configure the function settings of this pipeline segment
template<typename nodeType>
bool distMatrixPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
	/**
	    configPipe(std::map<std::string, std::string> &configMap)
	 
		@brief Configures the pipe and sets arguments based on the configMap passed. Called before execution (runPipe). If required values not found or configuration is invalid, returns false. 
		@tparam nodeType The data type of the simplex node.
		@param configMap The configuration map for this pipeline
        @return boolean
	*/
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();
		
	pipe = configMap.find("beta");
	if(pipe != configMap.end())
		this->beta = std::atof(configMap["beta"].c_str());
	
	pipe = configMap.find("betaMode");
	if(pipe != configMap.end())
		this->betaMode = configMap["betaMode"].c_str();
		
	this->ut = utils(strDebug, this->outputFile);
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		this->enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else return false;

	this->configured = true;
	this->ut.writeDebug("distMatrixPipe","Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");
	
	return true;
}

template<typename nodeType>
void distMatrixPipe<nodeType>::outputData(pipePacket<nodeType> &inData){
    /**
	    outputData(pipePacket<nodeType> &inData)
	 
		@brief Outputs distance matrix to a file if debug mode is true. 
		@tparam nodeType The data type of the simplex node.
		@param inData The pipePacket data being used in the pipeline.
	*/
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
template class distMatrixPipe<witnessNode>;
