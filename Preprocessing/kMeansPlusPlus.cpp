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
#include "kMeansPlusPlus.hpp"
#include "utils.hpp"

// basePipe constructor
kMeansPlusPlus::kMeansPlusPlus(){
	procName = "k-means++";
    return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket kMeansPlusPlus::runPreprocessor(pipePacket inData){
    //Arguments - num_clusters, num_iterations

    utils ut;
    std::vector<std::vector<double>> centroids;     //Storing centroids
    std::vector<int> labels;                        //Storing labels for mapping data to centroids

    //Initialize centroids (Plus plus mechanism with kmeans - Hartigan, Wong)
    
    //initialize first random centroid from data
    
    //mersenne twister random algorithm - used to have reproducible results from seed
    //  This seed should be recorded to reproduce after a run
    //  There may be multiple seeds in a run depending on how many times k-means is used
    static std::random_device seed;  
    static std::mt19937 gen(seed()); 
    
    std::uniform_int_distribution<size_t> distribution(0, inData.workData.originalData.size()-1);
    int index = distribution(gen);

    
    std::vector<std::vector<double>> tempCentroids; 
    std::vector<double> center_initial = inData.workData.originalData[index];
    centroids.push_back(center_initial);  //adding first mean to centroids
    tempCentroids.push_back(center_initial);
    std::vector<double> tempDist;



	//Determining initial centroids by probability / distance from previous centroid
	// 	Do we need to compare the picked centroid to all centroids previously picked 
	//		-> in order to maximize the distance/difference between centroids
     for(unsigned k = 0; k<num_clusters-1; k++){ //adding means 2 -> k-1 to centroids based on distance 
        int index_next = distribution(gen);
        bool picked = false;
        std::vector<double>  center_next = inData.workData.originalData[index_next];
            for(unsigned j=0; j<inData.workData.originalData.size()-1; j++) {
                auto dist_next = ut.vectors_distance(inData.workData.originalData[j], center_next);
                auto centers_dist = ut.vectors_distance(center_initial, center_next);
                auto probability = (centers_dist*centers_dist)/dist_next;
                
                bool valid_center = (rand() %100) < (probability) * 10;
                
                if(valid_center && (std::find(tempCentroids.begin(), tempCentroids.end(), inData.workData.originalData[index_next]) == tempCentroids.end())){
                    tempCentroids.push_back(inData.workData.originalData[index_next]);
                    picked = true;
                    break;
                }
            
            }
    }       
    
    int i = 0;
    while(tempCentroids.size() < num_clusters){
		if(std::find(tempCentroids.begin(), tempCentroids.end(), inData.workData.originalData[i]) == tempCentroids.end()){
			tempCentroids.push_back(inData.workData.originalData[i]);
		}
		i++;
	}
    
    
    //Now, we need to iterate over the centroids and classify each point
    //	Compute the new centroid as the geometric mean of the classified points
    //		IF the centroid has 0 points classified to it, we need to repick/reapproach (find out how)
    //	Replace old centroids with new centroids, rinse, repeat until convergence
    //		Convergence on # iterations, or a minimum movement of centroids (< .01)
    // OUTPUT: Final centroids, labeled data, sum of vectors in a label, count of points in label
    
    //Need to store: 
    //	Counts, the number of clusters in each classification
    //	summedClusters, the total cluster distance (for r_max and r_avg)
    //	summedCentroidVectors, for the geometric sum of the centroid
    //	lastCentroids, to track the change in cluster WCSSE
    std::vector<unsigned> lastLabels;
    
    //Iterate until we reach max iderations or no change in the classification
	for (int z = 0; z < num_iterations; z++){
		
		std::vector<double> summedClusters(num_clusters, 0); 
		std::vector<std::vector<double>> summedCentroidVectors(num_clusters, std::vector<double>(inData.workData.originalData[0].size(), 0));
		std::vector<double> counts(num_clusters, 0);
		std::vector<unsigned> curLabels;
		
		//For each point, classify it to a centroid
		for (unsigned j = 0; j < inData.workData.originalData.size() - 1; j++){
			double minDist = std::numeric_limits<double>::max();
			unsigned clusterIndex = 0;
			
			//Check each centroid for the minimum distance
			for (unsigned c = 0; c < tempCentroids.size(); c++){
				auto curDist = ut.vectors_distance(inData.workData.originalData[j], tempCentroids[c]);
				
				
				if(curDist < minDist){
					clusterIndex = c;
					minDist = curDist;
				}
			}
			
			for(int d = 0; d < inData.workData.originalData[j].size(); d++){
				summedCentroidVectors[clusterIndex][d] += inData.workData.originalData[j][d];
			}
			curLabels.push_back(clusterIndex);
			counts[clusterIndex] ++;
			summedClusters[clusterIndex] += minDist;
			
		}		
		
		//Otherwise, 
		//		Shift the centroid geometric centers to new centroids
		for(int i = 0; i < summedCentroidVectors.size(); i++){
			for(int d =0; d < summedCentroidVectors[0].size(); d++){
				summedCentroidVectors[i][d] = summedCentroidVectors[i][d] / counts[i];
			}
		}
		
		tempCentroids = summedCentroidVectors;
		
		//Check for convergence
		if(curLabels == lastLabels){
			break;
		}		
		
		
		lastLabels = curLabels;
		
	}
    
	std::cout << "Clustered data..." << std::endl;
	//Assign to the pipepacket
	inData.workData.originalData = tempCentroids;
	inData.workData.originalLabels = lastLabels;
	return inData;
}

// configPipe -> configure the function settings of this pipeline segment
bool kMeansPlusPlus::configPreprocessor(std::map<std::string, std::string> configMap){
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
