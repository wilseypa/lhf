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
        std::vector<double>  center_next = inData.workData.originalData[index_next];
            for(unsigned j=0; j<inData.workData.originalData.size()-1; j++) {
                auto dist_next = ut.vectors_distance(inData.workData.originalData[j], center_next);
                auto centers_dist = ut.vectors_distance(center_initial, center_next);
                auto probability = (centers_dist*centers_dist)/dist_next;
                
                bool valid_center = (rand() %100) < (probability) * 10;
                
                if(valid_center && (std::find(tempCentroids.begin(), tempCentroids.end(), inData.workData.originalData[index_next]) == tempCentroids.end())){
                    tempCentroids.push_back(inData.workData.originalData[index_next]);
                }
            
            }
    }       
    
    //Now, we need to iterate over the centroids and classify each point
    //	Compute the new centroid as the geometric mean of the classified points
    //		IF the centroid has 0 points classified to it, we need to repick/reapproach (find out how)
    //	Replace old centroids with new centroids, rinse, repeat until convergence
    //		Convergence on # iterations, or a minimum movement of centroids (< .01)
    // OUTPUT: Final centroids, labeled data, sum of vectors in a label, count of points in label
    
    std::vector<double> counts(num_clusters, 0);
    std::vector<double> summed_clusters(num_clusters, 0); 
    
	for (int z = 0; z < num_iterations; z++){
		int minDistanceLabel = 0;
		double minDist = std::numeric_limits<double>::max();
		for (int c = 0; c < tempCentroids.size(); c++){
			int minDistanceLabel = 0;
			
			for ( int j = 0; j < inData.workData.originalData.size(); j++){
				auto curDist = ut.vectors_distance(inData.workData.originalData[j], tempCentroids[c]);
				
				
				if(curDist < minDist){
					minDistanceLabel = c;
					minDist = curDist;
				}
				
				
			}
			
			counts[minDistanceLabel] ++;
			summed_clusters[c] += minDist;
			
			//Compute feature sum of distances for a label ** Nick to check on this and revise** 
			// summed_features[label] += ut.featureDistance(originalData, tempCentroids)
			
		}		
		
		//if z=num_iterations or converge: 
		//	Assign Final Labels
		//	Assign final centroid points (geometric center of classified points)
		//  Rmax, Ravg
		//	BREAK
		
		
		//Otherwise, 
		//		Shift the centroid geometric centers to new centroids
		
	}
   
   

	/*std::cout << tempCentroids.size() << " : " << tempCentroids[0].size() << std::endl;

	auto end = tempCentroids.end();  //removing duplicates from temp matrix
	for (auto it = tempCentroids.begin(); it != end; it++){
		end = std::remove(it +1, end, *it);
	}
	tempCentroids.erase(end, tempCentroids.end());
	
	
	std::cout << tempCentroids.size() << " : " << tempCentroids[0].size() << std::endl;
*/

	for(int i = 0; i<tempCentroids.size(); i++){
		centroids.push_back(tempCentroids[i]);  //adding unique centroids to centroids matrix
	} 


	for(unsigned i = 0; i < centroids.size(); i++){
	for(unsigned j = 0; j < centroids[j].size(); j++){
	std::cout << i + 1 << "->";
		std::cout << centroids[i][j]<< '\t';

	}
	std::cout << std::endl;

	}
	//labeling points so they can be assigned to a cluster
    //  This will eventually occur in parallel, should probably be removed (eventually)
    //  i.e. labels are not necessary to continue the PH computation but are for upscaling

    for(auto r = 0; r<num_iterations; r++){   //reassign points to centroids based on cmd arg iterations
        for(size_t point = 0; point<inData.workData.originalData.size(); point++){
            double dist_best = std::numeric_limits<double>::max(); //setting max starting distance
            //size_t cluster_best = 0;
            for(size_t cluster = 0; cluster<num_clusters; cluster++){
                const double cluster_dist = ut.vectors_distance(inData.workData.originalData[point], centroids[cluster]);
                if(cluster_dist < dist_best){
                    dist_best = cluster_dist; //assigning new best distance based on centroids
                    //cluster_best = centroids[cluster];
                    labels.push_back(cluster);  //labeling each point based on best cluster
                }   
            }
    
        }
    }


  /*for(unsigned i = 0; i < labels.size(); i++){
            std::cout << labels[i] << ",";
    }
    std::cout << "\n";*/  

    
    //assigning points to a cluster based on their label
    std::vector<std::vector<double>> assignments;
   // std::vector<double> assignments;
    int *ptr;
    for(int k=0; k<num_clusters; k++){  
        for(size_t point = 0; point<inData.workData.originalData.size(); point++){
            if(labels[point] == k){ 
             // std::cout << inData.workData.originalData[point][0] <<"point";
                assignments.push_back(inData.workData.originalData[point]);
            }
        }
        
    }  


 /*  for(unsigned i = 0; i < assignments.size(); i++){
		for(unsigned j = 0; j < assignments[i].size(); j++){
			std::cout << assignments[i][j] << '\t';
		}
		std::cout << std::endl;
	}   */

    //counting points in each cluster so mean can be computed
  
    for(size_t point=0; point<inData.workData.originalData.size(); point++){   
        summed_clusters.push_back(0);
        for(size_t j=0; j<num_clusters; j++){   
            if(labels[point] == j){
                counts[j] += 1;
                summed_clusters[point] += ut.vectors_distance(inData.workData.originalData[point], centroids[labels[point]]);
                
            }   
        }
    } 
	std::cout << std::endl;
	for(unsigned j = 0; j < counts.size(); j++){
		std::cout << counts[j] << '\t';
        std::cout << summed_clusters[j] << "\trAvg = " << summed_clusters[j]/counts[j] << std::endl;
	} 


    //summing points in each cluster so mean can be computed... not working yet

    
    std::cout << "Clustered data..." << std::endl;
    //Assign to the pipepacket
    inData.workData.originalData = centroids;
    inData.workData.originalLabels = labels;
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
