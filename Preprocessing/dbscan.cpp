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

// basePipe constructor
dbscan::dbscan(){
	procName = "dbscan";
    return;
}
//taking in preprocessor type

std::vector<int> dbscan::cluster(const std::vector<std::vector<double>> &data){
    std::vector<int> labels(data.size(), 0);
    kdTree tree(data);
    int clusterLabel = 0;

    int i=0;
    for(int i=0; i<data.size(); i++){
        if(labels[i] != 0) continue;

        std::vector<size_t> neighbors = tree.neighborhoodIndices(data[i], epsilon);
        if(neighbors.size() < minPoints) labels[i] = -1;
        else{
            clusterLabel++;
            expandCluster(data, labels, neighbors, clusterLabel, tree);
        }
    }

    return labels;
}

void dbscan::expandCluster(const std::vector<std::vector<double>> &data, 
                           std::vector<int> &labels, 
                           std::vector<size_t> &neighbors, 
                           int clusterLabel,
                           kdTree &tree){
    int i = 0;
    while(i < neighbors.size()){
        int pt = neighbors[i];
        if(labels[pt] == -1) labels[pt] = clusterLabel;
        else if(labels[pt] == 0){
            labels[pt] = clusterLabel;
            std::vector<size_t> ptNeighbors = tree.neighborhoodIndices(data[pt], epsilon);
            
            if(ptNeighbors.size() >= minPoints){
                neighbors.insert(neighbors.end(), ptNeighbors.begin(), ptNeighbors.end());
            }
        }
        ++i;
    }
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket dbscan::runPreprocessor(pipePacket inData){ //standalone preprocessor
    // inData.originalLabels = cluster(inData.originalData);
    return inData;
}

// configPipe -> configure the function settings of this pipeline segment
bool dbscan::configPreprocessor(std::map<std::string, std::string> configMap){
    std::string strDebug;
    
    auto pipe = configMap.find("debug");
    if(pipe != configMap.end()){
        debug = std::atoi(configMap["debug"].c_str());
        strDebug = configMap["debug"];
    }
    pipe = configMap.find("outputFile");
    if(pipe != configMap.end())
        outputFile = configMap["outputFile"].c_str();
    
    ut = utils(strDebug, outputFile);
    
    pipe = configMap.find("minPoints");
    if(pipe !=configMap.end())
        minPoints = std::stoi(configMap["minPoints"].c_str());
    else return false;

    pipe = configMap.find("epsilon");
    if(pipe != configMap.end())
        epsilon = std::stod(configMap["epsilon"].c_str());
    else return false;  
    
    ut.writeDebug("dbscan","Configured with parameters { minPoints: " + configMap["minPoints"] + ", epsilon: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");
}
