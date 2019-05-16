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
#include "coresetUtils.hpp"
#include "streamingUtils.hpp"
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


	int m = 5;
	int tracker = 0;
	 //open file stream  --> use file from input arg
	 //while open ... do
	std::vector<std::vector<double>> tempBucket;  //bucket 0
	std::vector<std::vector<double>> firstBucket; 
	std::vector<std::vector<std::vector<double>>> bucketManager; // holds all filled buckets
	 for(unsigned i = 0; i<inData.workData.originalData.size()-1; i++){
		
    
    
		for(unsigned b = 0; i<m-1; b++){
       tempBucket.push_back(inData.workData.originalData[i]);  //add first point to Bucket 1  
			 if(i%m == 0){
				 for(unsigned j = 0; j<tempBucket.size()-1; j++) {
					 bucketManager.push_back(tempBucket);   //filling 0th and 1st bucket
				 }
			 }                               
		}
		// if 0th and 1st bucket full, compute new coreset of size m from union of the points in b0 and b1
		//and put in b2, for b(i) = 2 -> numBuckets. then clear b0 so new points can be assigned
    







		}
    //for a data stream of n points, there are log2(n/m) + 2 buckets
		//extract coresets of size m from data stream -> merge and reduce technique - Har-Peled & Mazumdar 2004
		//new points inserted into Bucket B0 (0->m points)

       // choose next points "according to D^2" -> (dist(point, set)/(sum of[dist(set, point)])
		//when B0 full (now contains a coreset), move points to B1
    
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
