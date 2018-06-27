/*
 * distMatrix hpp + cpp extend the basePipe class for calculating the 
 * distance matrix from data input
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include "distMatrixPipe.hpp"
#include "utils.cpp"

// basePipe constructor
distMatrixPipe::distMatrixPipe(){

	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket distMatrixPipe::runPipe(pipePacket inData){
	utils ut;
	//Store our distance matrix
	std::vector<std::vector<double>> distMatrix;
	 
	//Iterate through each vector
	for(unsigned i = 0; i < inData.workData.workingData.size()-1; i++){
		//Grab a second vector to compare to 
		std::vector<double> temp;
		for(unsigned j = 0; j < inData.workData.workingData.size()-1; j++){

				//Calculate vector distance 
				auto dist = ut.vectors_distance(inData.workData.workingData[i],inData.workData.workingData[j]);
	
				temp.push_back(dist);	
		}
		distMatrix.push_back(temp);
	}
	
	//Assign to the pipePacket
	inData.workData.workingData = distMatrix;
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool distMatrixPipe::configPipe(std::map<std::string, std::string> configMap){
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	else return false;
	
	return true;
}

