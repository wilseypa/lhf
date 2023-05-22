/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute
 *
 */
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <random>
#include <chrono>
#include <string>
#include <numeric>
#include <iostream>
#include <functional> 
#include <vector>
#include "kMeansPlusPlus.hpp"
#include "cluster.hpp"
#include "utils.hpp"

/**
	@class kMeansPlusPlus
	@brief A class for implementing the k-means++ algorithm
	@tparam nodeType The type of data points to cluster
*/


/**
@brief kMeansPlusPlus class constructor

Initializes the kMeansPlusPlus class instance.
*/
template<typename nodeType>
kMeansPlusPlus<nodeType>::kMeansPlusPlus(){
	this->procName = "k-means++";
    return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment

/**
@brief Runs the configured functions of this pipeline segment
@param inData Input data to be processed by the preprocessor
*/
template<typename nodeType>
void kMeansPlusPlus<nodeType>::runPreprocessor(pipePacket<nodeType>& inData){
	if(!this->configured){
		this->ut.writeLog(this->procName,"Preprocessor not configured");
		return;
	}
	
    auto algo = kmeansplusplus();
    algo.clusterData(inData.workData, inData.workData, inData.centroidLabels, this->num_clusters, this->num_iterations, this->seed);
    
    
	return;
}

// configPipe -> configure the function settings of this pipeline segment
/**
Cluster the input data using k-means++ algorithm.

@param inData The input data to be clustered.

@return The clustered data in a pipePacket.
*/
template<typename nodeType>
bool kMeansPlusPlus<nodeType>::configPreprocessor(std::map<std::string, std::string>& configMap){
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = (std::atoi(configMap["debug"].c_str()) > 0 ? true : false);
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();
	
	this->ut = utils(strDebug, this->outputFile);
    pipe = configMap.find("clusters");
    if(pipe !=configMap.end())
        this->num_clusters = std::atoi(configMap["clusters"].c_str());
    else return false;
    
    pipe = configMap.find("seed");
    if(pipe != configMap.end())
		this->seed = std::atoi(configMap["seed"].c_str());

    pipe = configMap.find("iterations");
	if(pipe != configMap.end())
		this->num_iterations = std::atoi(configMap["iterations"].c_str());
	else return false;
	
	this->configured = true;
	this->ut.writeDebug(this->procName,"Configured with parameters { clusters: " + configMap["clusters"] + ", iterations: " + configMap["iterations"] + ", debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

//Explicit Template Class Instantiation
template class kMeansPlusPlus<simplexNode>;
template class kMeansPlusPlus<alphaNode>;
template class kMeansPlusPlus<witnessNode>;
