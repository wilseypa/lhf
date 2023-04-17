/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute
 *
 */

/**
 * @file streamingUtils.hpp
 * @brief Header file for the StreamingUtils class.
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
#include "streamingUtils.hpp"
#include "utils.hpp"
////// use these as post processing step for StreamKM on very large data sets 
/**
 * @brief The streamingUtils class is a utility class designed for use with StreamKM algorithm.
 * It contains post-processing steps for StreamKM on very large data sets.
 *
 * @tparam nodeType The data type of the elements stored in the streamingUtils class.
 */

/**
 * @brief Default constructor for the streamingUtils class.
 * Initializes an empty streamingUtils object.
 *
 * @tparam nodeType The data type of the elements stored in the streamingUtils class.
 */

// basePipe constructor
template <typename nodeType>
streamingUtils<nodeType>::streamingUtils(){
  
}

/**
 * @brief Performs k-means clustering on a dataset.
 *
 * This class member function takes in a dataset and performs k-means clustering
 * using a specified number of clusters and iterations.
 *
 * @tparam nodeType Template parameter for the data type of elements in the dataset.
 * @param kHat Input dataset, a vector of vector of doubles representing data points.
 * @return A vector of vector of doubles representing the centroids of the clusters.
 */

template <typename nodeType>
std::vector<std::vector<double>> streamingUtils<nodeType>::kMeans(std::vector<std::vector<double>>& kHat){ //"batch"/ normal kmeans clustering to be performed on data
    utils ut;
     int numIterations;
     int numClusters;
    
    static std::random_device seed;  
    static std::mt19937 gen(seed()); 
    std::vector<std::vector<double>> tempCentroids; 
    std::vector<std::vector<double>> tempClusters; 
    std::uniform_int_distribution<size_t> distribution(0, kHat[0].size());
     std::vector<unsigned> lastLabels;

    for(int i = 0; i<numClusters; i++){   //randomly picking centroids from kHat
        int index = distribution(gen);
        std::vector<double> centroid = kHat[index];
        tempCentroids.push_back(centroid);  //adding random means  to tempCentroids
    }
    


    for(int z = 0; z < numIterations; z++){    //iterate until max iteration or no change in classification
        std::vector<double> summedClusters(numClusters, 0);
        std::vector<std::vector<double>> summedCentroidVectors(numClusters, std::vector<double>(kHat[0].size(), 0));
		std::vector<double> counts(numClusters, 0);
		std::vector<unsigned> curLabels;


       //classify each point to a centroid
        for(int j = 0; j< kHat.size()-1; j++){
            double minDist =  std::numeric_limits<double>::max();
            unsigned clusterIndex = 0;

            //check each centroid for the minimum distance
            for(int c = 0; c<tempCentroids.size(); c++){
                auto curDist = ut.vectors_distance(kHat[j], tempCentroids[c]);

                    if(curDist < minDist){
                        clusterIndex = c;
                        minDist = curDist;
                    }
            }

            for(int d = 0; d <  kHat[j].size(); d++){
				summedCentroidVectors[clusterIndex][d] += kHat[j][d];
			}
			curLabels.push_back(clusterIndex);
			counts[clusterIndex] ++;
			summedClusters[clusterIndex] += minDist;
        }

        tempClusters = summedCentroidVectors;

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
    
    return tempClusters;

}


/**
 * @brief Performs ball k-means clustering on a dataset.
 *
 * This class member function takes in a dataset and performs ball k-means clustering
 * using a specified number of clusters, radius and centroids.
 *
 * @tparam nodeType Template parameter for the data type of elements in the dataset.
 * @param clusters Input dataset, a vector of vector of doubles representing data points.
 * @return A vector of vector of doubles representing the centroids of the clusters.
 */

template <typename nodeType>
std::vector<std::vector<double>> streamingUtils<nodeType>::ballKmeans(std::vector<std::vector<double>>& clusters){

								   //select points in each cluster closest to each centroid within ball radius
									 //compute center of mass of those points --> becomes final centroid 
    utils ut;
    static std::random_device seed;  
    static std::mt19937 gen(seed()); 
    std::vector<std::vector<double>> tempCentroids;
    std::uniform_int_distribution<size_t> distribution(0, clusters.size());
        for(int i = 0; i<2; i++){   //randomly picking  2 centroids from clusters
            int index = distribution(gen);
            std::vector<double> centroid = clusters[index];
            tempCentroids.push_back(centroid);  //adding random means  to tempCentroids
         }

    double ballRadius = ut.vectors_distance(tempCentroids[0], tempCentroids[1])/3;  //radius = (distance between 2 init centroids)/3
    std::vector<double> radius;
    radius.push_back(ballRadius);
    std::vector<std::vector<double>> summedCentroidVectors;
    std::vector<double> counts;
    std::vector<unsigned> curLabels;
    std::vector<double> summedClusters;
    std::vector<unsigned> lastLabels;

    for(int i =0; i< clusters.size()-1; i++){
        unsigned clusterIndex = 0;
        //checking each point if it lies within ball radius
        	for (unsigned c = 0; c < tempCentroids.size(); c++){
				auto curDist = ut.vectors_distance(clusters[i], radius);
				
				if(curDist < ballRadius){
					clusterIndex = c;
				}
			}

            for(int j = 0; j<clusters.size(); j++){
                	summedCentroidVectors[clusterIndex][j] += clusters[i][j];
			}
			curLabels.push_back(clusterIndex);
			counts[clusterIndex] ++;
	

    }


    for(int i = 0; i < summedCentroidVectors.size(); i++){  //moving each centroid to its center of mass per cluster
		for(int d =0; d < summedCentroidVectors[0].size(); d++){
			summedCentroidVectors[i][d] = summedCentroidVectors[i][d] / counts[i];
		}
	}
		
		tempCentroids = summedCentroidVectors;
		
		//Check for convergence
	/*	if(curLabels == lastLabels){
			break;
		}	*/

    return tempCentroids;
    }


   
