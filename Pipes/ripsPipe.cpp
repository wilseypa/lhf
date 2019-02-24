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
	pipeType = "Rips_Inductive";
	return;
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket ripsPipe::runPipe(pipePacket inData){
	utils ut;
	
	inData.workData.complex->expandDimensions(dim);
		
	std::cout << "\tComplex Size: " << inData.workData.complex->simplexCount() << std::endl;
	std::cout << "\tComplex Mem: " << inData.workData.complex->getSize() << std::endl;
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool ripsPipe::configPipe(std::map<std::string, std::string> configMap){
	
	
	auto pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	
	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void ripsPipe::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");
	
	for (int i = 0; i < inData.workData.complex->weightedGraph.size(); i++){
		for(auto a : inData.workData.complex->weightedGraph[i]){
			for(auto d : a){
				file << d << ",";
			}
			file << "\n";
		}
	}
	
	file.close();
	return;
}

