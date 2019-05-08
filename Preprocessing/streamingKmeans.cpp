/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute
 *
 */
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <random>
#include <chrono>
#include <string>
#include <numeric>
#include <iostream>
#include <functional> 
#include <vector>
#include "streamingKmeans.hpp"
#include "utils.hpp"
// basePipe constructor
streamingKmeans::streamingKmeans(){
	procName = "streamingKmeans";
	return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket streamingKmeans::runPreprocessor(pipePacket inData){
	//Arguments - num_clusters, num_iterations

	utils ut;
	static std::random_device seed;
	static std::mt19937 gen(seed());
	std::uniform_int_distribution<size_t> distribution(0, inData.workData.originalData.size()-1);
	 //open file stream  --> use file from input arg
	 //while open ... do
		//extract coresets of size m from data stream -> merge and reduce technique - Har-Peled & Mazumdar 2004

		//new points inserted into Bucket B0 (0->m points)

		//when B0 full, move points to B1
			//if B1 full (contains m points) make new coreset from B0 and B1 and move to B2 .. repeat for i







	 //close file stream

	// do kmeans++ on reduced dataset (reduce step)

	 //return clustered data 
	
	
	
	
     



	return inData;
}

// configPipe -> configure the function settings of this pipeline segment
bool streamingKmeans::configPreprocessor(std::map<std::string, std::string> configMap){
  auto preprocessor = configMap.find("clusters");
    if(preprocessor !=configMap.end())
        num_clusters = std::atoi(configMap["clusters"].c_str());
    else return false;

    preprocessor = configMap.find("iterations");
	if(preprocessor != configMap.end())
		num_iterations = std::atoi(configMap["iterations"].c_str());
	else return false;	
	
	
	return true;
}

//assigning points to nearest cluster
