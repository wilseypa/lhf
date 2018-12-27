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
#include "kMeansPlusPlus.hpp"
#include "utils.hpp"
// basePipe constructor
kMeansPlusPlus::kMeansPlusPlus(){
	return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket kMeansPlusPlus::runPreprocessor(pipePacket inData){
	//Arguments - How many clusters? Convergence value
				utils ut;
				//Storing centroids
				std::vector<std::vector<double>> centroids;
				//Storing labels for mapping data to centroids
				std::vector<int> labels;



	//Initialize centroids (Plus plus mechanism with kmeans - Hartigan, Wong)
    
	//initialize first random centroid from data
	static std::random_device seed;  
	static std::mt19937 gen(seed()); //mersenne twister random algorithm
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

	int num_clusters = 5; //need to set out of cmd line later... testing for now

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
	//assigning points to centroids
	int num_iterations = 250; //pass in from cmd line? using this to test for now
	std::vector<size_t> assignments(inData.workData.originalData.size());
	for (size_t point = 0; point<inData.workData.originalData.size(); point++){
		double dist_best = std::numeric_limits<double>::max(); //setting max starting distance
		size_t cluster_best = 0;
		for(size_t cluster = 0; cluster <num_clusters; cluster++){
			const double cluster_dist = ut.vectors_distance(inData.workData.originalData[point], centroids[cluster]);
            if(cluster_dist < dist_best){
				dist_best = cluster_dist; //assigning new best distance based on centroids
				//cluster_best = centroids[cluster];
				labels.push_back(cluster);
			}
		}
	
		assignments[point] = cluster_best;
	}
    
	
	std::cout << "No run function defined for: " << procName << std::endl;
	//Assign to the pipepacket
    inData.workData.originalData = centroids;

	return inData;
}

// configPipe -> configure the function settings of this pipeline segment
bool kMeansPlusPlus::configPreprocessor(std::map<std::string, std::string> configMap){

	std::cout << "No configure function defined for: " << procName << std::endl;
    
	return true;
}

//assigning points to nearest cluster
