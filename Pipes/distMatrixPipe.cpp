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
#include <algorithm>
#include "distMatrixPipe.hpp"
#include "utils.hpp"

// basePipe constructor
distMatrixPipe::distMatrixPipe(){
	pipeType = "DistMatrix";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket distMatrixPipe::runPipe(pipePacket inData){
	utils ut;
	//Store our distance matrix
	std::vector<std::vector<double>> distMatrix;
	 
	//Iterate through each vector
	for(unsigned i = 0; i < inData.workData.originalData.size()-1; i++){
		//Grab a second vector to compare to 
		std::vector<double> temp;
		for(unsigned j = 0; j < inData.workData.originalData.size()-1; j++){

				//Calculate vector distance 
				auto dist = ut.vectors_distance(inData.workData.originalData[i],inData.workData.originalData[j]);
	
				temp.push_back(dist);	
		}
		distMatrix.push_back(temp);
	}
	
	inData.workData.complex->setDistanceMatrix(distMatrix);
	
	std::cout << "\tDist Matrix Size: " << distMatrix.size() << " x " << distMatrix.size() << std::endl;
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

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void distMatrixPipe::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");
	
	for(auto a : inData.workData.complex->distMatrix){
		for(auto d : a){
			file << d << ",";
		}
		file << "\n";
	}
	
	file.close();
	return;
}

