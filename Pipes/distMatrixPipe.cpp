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
distMatrixPipe::distMatrixPipe(){
	pipeType = "DistMatrix";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket distMatrixPipe::runPipe(pipePacket inData){
	utils ut;
	
	//Store our distance matrix
	std::vector<std::vector<double>> distMatrix (inData.originalData.size(), std::vector<double>(inData.originalData.size(),0));
	
	//Iterate through each vector
	for(unsigned i = 0; i < inData.originalData.size(); i++){
		if(!inData.originalData[i].empty()){
		
			//Grab a second vector to compare to 
			std::vector<double> temp;
			for(unsigned j = i+1; j < inData.originalData.size(); j++){

					//Calculate vector distance 
					auto dist = ut.vectors_distance(inData.originalData[i],inData.originalData[j]);
					
					if(dist < maxEpsilon)
						inData.weights.insert(dist);
					distMatrix[i][j] = dist;
					//temp.push_back(dist);	
			}
			//distMatrix.push_back(temp);
		}
	}
	
	inData.complex->setDistanceMatrix(distMatrix);
	
	inData.weights.insert(0.0);
	inData.weights.insert(maxEpsilon);
	//std::sort(inData.weights.begin(), inData.weights.end(), std::greater<>());
	
	std::cout << "\tDist Matrix Size: " << distMatrix.size() << " x " << distMatrix.size() << std::endl;
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool distMatrixPipe::configPipe(std::map<std::string, std::string> configMap){
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	else return false;
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void distMatrixPipe::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");
	
	for(auto a : inData.complex->distMatrix){
		for(auto d : a){
			file << d << ",";
		}
		file << "\n";
	}
	
	file.close();
	return;
}

