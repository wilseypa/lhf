/*
 * betaSkeletonBasedComplexPipe hpp + cpp extend the basePipe class for calculating the 
 * beta Skeleton Based Complex generation for data input
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
#include "betaSkeletonBasedComplex.hpp"
#include "utils.hpp"

// basePipe constructor
betaSkeletonBasedComplexPipe::betaSkeletonBasedComplexPipe(){
	pipeType = "betaSkeletonBasedComplex";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
void betaSkeletonBasedComplexPipe::runPipe(pipePacket &inData){
	// Generate Beta Skeleton Based Complex
	ut.writeDebug("betaSkeletonBasedComplex Pipe", "\tbetaSkeletonBasedComplex Size: ");
	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool betaSkeletonBasedComplexPipe::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = configMap["outputFile"].c_str();
	
	pipe = configMap.find("beta");
	if(pipe != configMap.end())
		beta = std::atof(configMap["beta"].c_str());
		
	ut = utils(strDebug, outputFile);
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else return false;

	configured = true;
	ut.writeDebug("betaSkeletonBasedComplex Pipe ","Configured with parameters { eps: " + configMap["epsilon"] + configMap["beta"] + " , debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void betaSkeletonBasedComplexPipe::outputData(pipePacket &inData){
	// Output related to betaSkeletonBasedComplex
	return;
}

