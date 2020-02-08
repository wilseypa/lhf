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
#include "dbscan.hpp"
#include "utils.hpp"
#include "kdTree.hpp"
//////// DBSCAN algorithm for standalone clustering or as an initialization step for DenStream /////////

std::vector<int> dbscan::cluster(const std::vector<std::vector<double>> &data, int minPoints = 20, double epsilon = 16, int size = 0){
    if(size == 0){
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
        if(labels[pt] == -1) labels[pt] = clusterLabel;
        else if(labels[pt] == 0){ //Not yet visited
            labels[pt] = clusterLabel;
            std::vector<size_t> ptNeighbors = tree.neighborhoodIndices(data[pt], epsilon);
            
            if(ptNeighbors.size() >= minPoints){ //Also a core point
                neighbors.insert(neighbors.end(), ptNeighbors.begin(), ptNeighbors.end()); //Add all neighbors to cluster
            }
        }
        ++i;
    }
}