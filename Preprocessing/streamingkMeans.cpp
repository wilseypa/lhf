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
#include <vector>
#include "streamingkMeans.hpp"
#include "utils.hpp"
// basePipe constructor
streamingkMeans::streamingkMeans(){
	return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket streamingkMeans::runPreprocessor(pipePacket inData){
	//Arguments - num_clusters, num_iterations

	utils ut;
	std::vector<std::vector<double>> centroids;		//Storing centroids
	

    num_clusters = 5; //for testing purposes.... set out of cmd later
	//Initialize centroids (Plus plus mechanism with kmeans - Hartigan, Wong)
    
	//initialize first random centroid from data
	
	//mersenne twister random algorithm - used to have reproducible results from seed
	//	This seed should be recorded to reproduce after a run
	//	There may be multiple seeds in a run depending on how many times k-means is used
	static std::random_device seed;  
	static std::mt19937 gen(seed()); 
	
	std::uniform_int_distribution<size_t> distribution(0, inData.workData.originalData.size()-1);
    int index = distribution(gen);
	
	std::vector<double> center_initial = inData.workData.originalData[index];

	 //compute squared distances from initial center
	std::vector<double> tempdist;
    for(unsigned i = 0; i<inData.workData.originalData.size()-1; i++) {        
	    auto dist = ut.vectors_distance(inData.workData.originalData[i], center_initial);
		tempdist.push_back(dist);
	}
	
	//choose kth centroid based on probability of dist/(sum of dist)
    for(unsigned k = 0; k<num_clusters; k++){
        int index_next = distribution(gen);
		//std::vector<int> temp;
		std::vector<double>  center_next = inData.workData.originalData[index_next];
		    for(unsigned j=0; j<inData.workData.originalData.size()-1; j++) {
				auto dist_next = ut.vectors_distance(inData.workData.originalData[j], center_next);
				auto centers_dist = ut.vectors_distance(center_initial, center_next);
				auto probability = centers_dist/dist_next;
				
				bool valid_center = (rand() %100) < (probability) * 10;
				
				if(valid_center){
					centroids.push_back(inData.workData.originalData[index_next]);
				//	centroids.insert(std::end(centroids), std::begin(center_next), std::end(center_next));
					k++;
				}
				else{
				    k = k;
				}
			}
	} 


	std::vector<double> counts(num_clusters, 0);
	std::vector<double> new_mean;
	std::vector<std::vector<double>> old_mean;
	std::vector<std::vector<double>> temp_index;
	for (size_t i =0; i<inData.workData.originalData.size()-1; i++){
		for(size_t k = 0; k<num_clusters; k++){
			double dist_best = std::numeric_limits<double>::max(); //setting starting distance
			double tempdist = ut.vectors_distance(inData.workData.originalData[i], centroids[k]);
            if (tempdist <dist_best){  //finding the nearest centroid to current point
				dist_best = tempdist;  
				counts[k]++;
				old_mean.push_back(centroids[k]);
				temp_index.push_back(inData.workData.originalData[i]);
				new_mean[k] = (old_mean[k][0] + ( (1/counts[k])* sqrt(ut.vectors_distance(temp_index[k], old_mean[k])))); //modifying old mean as new points come in
				centroids.push_back(new_mean);
			}
		}


	}  
    
	
    std::cout << "Clustered data..." << std::endl;
	//Assign to the pipepacket
    inData.workData.originalData = centroids; 
     



	return inData;
}

// configPipe -> configure the function settings of this pipeline segment
bool streamingkMeans::configPreprocessor(std::map<std::string, std::string> configMap){

	num_clusters = stoi(configMap["clusters"]);			
	
	
	return true;
}

//assigning points to nearest cluster
