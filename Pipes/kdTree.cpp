/*
 * kd-Tree hpp + cpp extend the basePipe class for creating a k-d Tree for data input
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
#include "kdTree.hpp"
#include "utils.hpp"

// basePipe constructor
kdTreePipe::kdTreePipe(){
	pipeType = "kdTree";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
void kdTreePipe::runPipe(pipePacket &inData){
	// Create kd-Tree for input point cloud
	ut.writeDebug("kd-Tree", "\tkd-Tree Size: ");
	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool kdTreePipe::configPipe(std::map<std::string, std::string> &configMap){
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
	
	pipe = configMap.find("beta");
	if(pipe != configMap.end())
		beta = std::atof(configMap["beta"].c_str());
	else return false;

	configured = true;
	ut.writeDebug("kd-Tree Pipe","Configured with parameters { beta: " + configMap["beta"] + " , debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void kdTreePipe::outputData(pipePacket &inData){
	// kd-Tree Ouptput
	return;
}

