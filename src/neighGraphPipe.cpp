/*
 * distMatrix hpp + cpp extend the basePipe class for calculating the 
 * distance matrix from data input
 * 
 */

#include <string>
#include <iostream>
#include <iterator>
#include <vector>
#include <functional>
#include "neighGraphPipe.hpp"
#include "utils.hpp"


// basePipe constructor
neighGraphPipe::neighGraphPipe(){
	pipeType = "NeighborhoodGraph";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket neighGraphPipe::runPipe(pipePacket inData){	
	utils ut;
	//Store our nodes, edges, and weights
	std::vector<unsigned> nodeIndex(inData.workData.workingData.size()-1);
	std::vector<std::vector<unsigned>> edges;
	std::vector<double> weights;
	 
	 
	//Iterate through each vector
	for(unsigned i = 0; i < inData.workData.workingData.size(); i++){
		//Grab a second vector to compare to 
		for(unsigned j = 0; j < inData.workData.workingData.size()-i; j++){
			
			if (i != j+1){
				
				//Calculate vector distance 
				auto dist = ut.vectors_distance(inData.workData.workingData[i],inData.workData.workingData[j]);
				
				//Filter distances <= epsilon, > 0 (same point)
				if(dist <= epsilon && dist > 0){
					std::vector<unsigned> edge = {i,j};
					edges.push_back(edge);
					weights.push_back(dist);
				}
			}
		
		}
	}
	inData.workData.edges = edges;
	inData.workData.weights = weights;
	
	/*
	std::cout << "test\t" << inData.workData.workingData[0].size() << std::endl;
	std::cout << nodeIndex.size() << "\t" << edges.size() << "\t" << weights.size() << std::endl << std::endl << std::endl;
	std::cout << std::endl << std::endl;
	for(unsigned i = 0; i < edges.size(); i++){
		for (unsigned j = 0; j < edges[i].size(); j++){
			std::cout << edges[i][j] << "\t";
		}
		std::cout << weights[i] << std::endl;		
	}
	*/
	
	
	return inData;
}



// configPipe -> configure the function settings of this pipeline segment
bool neighGraphPipe::configPipe(std::map<std::string, std::string> configMap){
	auto pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		epsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	else return false;
	
	return true;
}

