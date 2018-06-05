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

// basePipe constructor
distMatrixPipe::distMatrixPipe(){

	return;
}


// Calculate the euclidean distance between two vectors
double vectors_distance(const std::vector<double>& a, const std::vector<double>& b){
		std::vector<double> temp;
		
		if(b.size() == 0)
			return 0;
		
		std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp),[](double e1, double e2) {return pow((e1-e2),2);});	
		return sqrt(std::accumulate(temp.begin(), temp.end(), 0.0));
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket distMatrixPipe::runPipe(pipePacket inData){
	
	//Store our distance matrix
	std::vector<std::vector<double>> distMatrix;
	 
	//Iterate through each vector
	for(unsigned i = 0; i < inData.workData.workingData.size()-1; i++){
		//Grab a second vector to compare to 
		std::vector<double> temp;
		for(unsigned j = 0; j < inData.workData.workingData.size()-1; j++){

				//Calculate vector distance 
				auto dist = vectors_distance(inData.workData.workingData[i],inData.workData.workingData[j]);
	
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
	
	return true;
}

