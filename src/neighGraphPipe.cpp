/*
 * distMatrix hpp + cpp extend the basePipe class for calculating the 
 * distance matrix from data input
 * 
 */

#include <string>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <functional>
#include "neighGraphPipe.hpp"


// basePipe constructor
neighGraphPipe::neighGraphPipe(){

	return;
}

double vectors_distance(const std::vector<double>& a, const std::vector<double>& b)
{
		
		std::vector<double> temp;
		
		if(b.size() == 0)
			return 0;
		
		std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp),[](double e1, double e2) {return pow((e1-e2),2);});
	
		return sqrt(std::accumulate(temp.begin(), temp.end(), 0.0));
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket neighGraphPipe::runPipe(pipePacket inData){	
	
	//Store our nodes, edges, and weights
	std::vector<int> nodeIndex(inData.workData.workingData.size()-1);
	std::vector<std::vector<int>> edges;
	std::vector<double> weights;
	 
	 
	//Iterate through each vector
	for(unsigned i = 0; i < inData.workData.workingData.size(); i++){
		//Grab a second vector to compare to 
		for(unsigned j = 0; j < inData.workData.workingData.size()-i; j++){
			
			if (i != j+1){
				
				//Calculate vector distance 
				auto dist = vectors_distance(inData.workData.workingData[i],inData.workData.workingData[j]);
	
				//Filter distances <= epsilon, > 0 (same point)
				if(dist <= epsilon && dist > 0){
					std::vector<int> edge = {i,j};
					edges.push_back(edge);
					weights.push_back(dist);
				}
			}
		
		}
	}
	
	//std::cout << nodeIndex.size() << "\t" << edges.size() << "\t" << weights.size() << std::endl;
	//td::cout << std::endl << std::endl;
	//for(unsigned i = 0; i < edges.size(); i++){
	//	for (unsigned j = 0; j < edges[1].size(); j++){
	//		std::cout << edges[i][j] << "\t";
	//	}
	//	std::cout << weights[i] << std::endl;		
	//}
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool neighGraphPipe::configPipe(std::map<std::string, std::string> configMap){
	
	auto pipe = configMap.find("epsilon");
	if(pipe != configMap.end()){
		epsilon = std::atof(configMap["epsilon"].c_str());
	}
	
	return true;
}

