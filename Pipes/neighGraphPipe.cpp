/*
 * distMatrix hpp + cpp extend the basePipe class for calculating the 
 * distance matrix from data input
 * 
 */

#include <string>
#include <iostream>
#include <fstream>
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

	//Iterate through each vector, inserting into simplex storage
	for(unsigned i = 0; i < inData.workData.originalData.size(); i++){
		
		//insert data into the complex (SimplexArrayList, SimplexTree)
		inData.workData.complex->insert(inData.workData.originalData[i]);
	}	


	std::cout << "returning..." << std::endl;
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


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void neighGraphPipe::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");
	
	for(auto a : inData.workData.complex->weightedGraph[1]){
		for(auto d : a){
			file << d << ",";
		}
		file << "\n";
	}
	
	file.close();
	return;
}
