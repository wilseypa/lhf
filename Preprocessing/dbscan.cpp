/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute
 *
 */

/**

    @file
    @brief This file includes several C++ standard libraries, as well as user-defined header files.
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
#include "dbscan.hpp"
#include "utils.hpp"
#include "kdTree.hpp"
//////// DBSCAN algorithm for standalone clustering or as an initialization step for DenStream /////////

/**

    @brief This function clusters the input data using DBSCAN algorithm.

    @param data The input data to cluster

    @param minPoints The minimum number of points required to form a dense region (default value: 20)

    @param epsilon The maximum distance between two points for them to be considered as neighbors (default value: 16)

    @param size The number of points to run DBSCAN on (default value: 0, which means the whole data set is clustered)

    @return A vector of cluster labels for each data point. Noise points are labeled as -1.
*/
std::vector<int> dbscan::cluster(const std::vector<std::vector<double>> &data, int minPoints = 20, double epsilon = 16, int size = 0){
    //Run DBSCAN on the first size points - to initiate DenStream
    if(size == 0){ //If not specified, cluster whole data set
        size = data.size();
    }

    std::vector<int> labels(size, 0);
    kdTree tree(data, size); //KDTree for efficient nearest neighbor search
    int clusterLabel = 0;

    for(int i=0; i<size; i++){
        if(labels[i] != 0) continue; //Already visited

        std::vector<size_t> neighbors = tree.neighborhoodIndices(data[i], epsilon); //All neighbors in epsilon-ball
        if(neighbors.size() < minPoints) labels[i] = -1; //Noise
        else{
            clusterLabel++;
            expandCluster(data, labels, neighbors, clusterLabel, tree, minPoints, epsilon); //Create new cluster
        }
    }

    return labels;
}

/**

    @brief Expands a cluster in the DBSCAN algorithm.
    @param data The data set to be clustered.
    @param labels The cluster labels of the data points.
    @param neighbors The indices of the neighbors of the current point.
    @param clusterLabel The label of the current cluster.
    @param tree The KDTree used for efficient nearest neighbor search.
    @param minPoints The minimum number of points required to form a cluster.
    @param epsilon The radius of the neighborhood used in the clustering.
*/

void dbscan::expandCluster(const std::vector<std::vector<double>> &data, 
                           std::vector<int> &labels, 
                           std::vector<size_t> &neighbors, 
                           int clusterLabel,
                           kdTree &tree,
                           int minPoints, 
                           double epsilon){
    int i = 0;
    while(i < neighbors.size()){
        int pt = neighbors[i];
        if(labels[pt] == -1) labels[pt] = clusterLabel; //Add noise points to cluster
        else if(labels[pt] == 0){ //Not yet visited
            labels[pt] = clusterLabel;
            std::vector<size_t> ptNeighbors = tree.neighborhoodIndices(data[pt], epsilon);
            
            if(ptNeighbors.size() >= minPoints){ //Also a core point
                neighbors.insert(neighbors.end(), ptNeighbors.begin(), ptNeighbors.end()); //Add all neighbors of core point to queue to traverse
            }
        }
        ++i;
    }
}