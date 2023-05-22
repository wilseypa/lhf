/**
 * @file utils.hpp
 *
 * @brief Contains utility functions for the LFH system (https://github.com/wilseypa/LFH).
 */


#include "cluster.hpp"
#include <cstdlib>
#include <chrono>
#include <random>



void kmeansplusplus::clusterData(std::vector<std::vector<double>> &orig, std::vector<std::vector<double>> &cent, std::vector<unsigned> &retLabels,  int num_clusters, int num_iterations, int seed = -1){
    //Arguments - num_clusters, num_iterations
    std::vector<int> labels; //Storing labels for mapping data to centroids
    unsigned dim = orig[0].size();

    //Initialize centroids (Plus plus mechanism with kmeans - Hartigan, Wong)
    
    //initialize first random centroid from data
    
    //mersenne twister random algorithm - used to have reproducible results from seed
    //  This seed should be recorded to reproduce after a run
    //  There may be multiple seeds in a run depending on how many times k-means is used
    //static std::random_device seed;
    if(seed < 0){
		seed = time(NULL);
    }
    static std::mt19937 gen(seed); 
    
    std::uniform_int_distribution<size_t> distribution(0, orig.size()-1);
    int index = distribution(gen); //Randomly choose the first centroid

    std::vector<std::vector<double>> centroids(num_clusters, std::vector<double>(dim, 0)); 
    std::vector<double> dist(orig.size(), std::numeric_limits<double>::max()); //Distance to nearest centroid

    centroids[0] = orig[index]; //Adding first mean to centroids

	//Determining initial centroids by probability / distance from previous centroid
	// 	Do we need to compare the picked centroid to all centroids previously picked 
	//		-> in order to maximize the distance/difference between centroids
    for(unsigned k = 1; k<num_clusters; k++){ //adding means 2 -> k to centroids based on distance 
     	double maxDist = 0; //Choose the point that is farthest from its closest centroid
     	int clusterIndex = 0;

        for(unsigned j=0; j<orig.size(); j++) {

			double curDist = this->ut.vectors_distance(orig[j], centroids[k-1]);	
			if(curDist < dist[j]) dist[j] = curDist;

			if(dist[j] > maxDist){
				maxDist = dist[j];
				clusterIndex = j;
			}
     	}

        centroids[k] = orig[clusterIndex];
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
	std::vector<std::vector<double>> summedCentroidVectors(num_clusters, std::vector<double>(dim, 0));

    
    //Iterate until we reach max iderations or no change in the classification
	for (int z = 0; z < num_iterations; z++){		
		std::vector<unsigned> counts(num_clusters, 0);
		std::vector<unsigned> curLabels(orig.size());
		
		//For each point, classify it to a centroid
		for (unsigned j = 0; j < orig.size(); j++){
			double minDist = std::numeric_limits<double>::max();
			unsigned clusterIndex = 0;
			
			//Check each centroid for the minimum distance
			for (unsigned c = 0; c < centroids.size(); c++){
				auto curDist = this->ut.vectors_distance(orig[j], centroids[c]);
				
				if(curDist < minDist){
					clusterIndex = c;
					minDist = curDist;
				}
			}
			
			if(z == 0){
				for(int d = 0; d < dim; d++)
					summedCentroidVectors[clusterIndex][d] += orig[j][d];
			} else if(lastLabels[j] != clusterIndex){
				for(int d = 0; d < dim; d++){
					summedCentroidVectors[lastLabels[j]][d] -= orig[j][d];
					summedCentroidVectors[clusterIndex][d] += orig[j][d];
				}
			}
			curLabels[j] = clusterIndex;
			counts[clusterIndex]++;			
		}		
		
		//Otherwise, 
		//		Shift the centroid geometric centers to new centroids
		for(int i = 0; i < summedCentroidVectors.size(); i++){
			if(counts[i] == 0){
				centroids[i] = orig[distribution(gen)];
			} else{
				for(int d = 0; d < summedCentroidVectors[0].size(); d++){
					centroids[i][d] = summedCentroidVectors[i][d] / counts[i];
				}
			}
		}
		
		//Check for convergence
		if(curLabels == lastLabels) break;
		lastLabels = curLabels;
		
	}
    
    
	//Assign to the pipepacket
	cent = centroids;
	retLabels = lastLabels;
    
    
    
}
