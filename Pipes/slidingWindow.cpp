/*
 * streamVR hpp + cpp extend the basePipe class for 
 * 
 */

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include <functional>
#include "slidingWindow.hpp"
#include "utils.hpp"


// basePipe constructor
slidingWindow::slidingWindow(){
	pipeType = "SlidingWindow";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket slidingWindow::runPipe(pipePacket inData){	
	utils ut;	
	
	//Iterate through each vector, inserting into simplex storage
	for(unsigned i = 0; i < inData.originalData.size(); i++){
		if(!inData.originalData[i].empty()){
			//insert data into the complex (SimplexArrayList, SimplexTree)
			inData.complex->insert(inData.originalData[i]);
			
		}
	}	
	
	
	inData.complex->expandDimensions(dim);

	return inData;
}



// configPipe -> configure the function settings of this pipeline segment
bool slidingWindow::configPipe(std::map<std::string, std::string> configMap){
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
	
	ut.writeDebug("slidingWindow","Configured with parameters { dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void slidingWindow::outputData(pipePacket inData){
	std::ofstream file ("output/" + pipeType + "_output.csv");
	
	if(inData.complex->simplexType == "simplexArrayList"){
		for(auto a : inData.complex->weightedGraph[1]){
			for(auto d : a){
				file << d << ",";
			}
			file << "\n";
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
