#pragma once

// Header file for naiveWindow class - see slidingWindow.cpp for descriptions
#include <map>
#include <unordered_map>
#include "basePipe.hpp"
#include "utils.hpp"

template <typename nodeType>
class slidingWindow : public basePipe<nodeType> {
private:
    double epsilon;
    int dim;
    int repCounter = 0;
    std::string inputFile;
    void runSubPipeline();
    std::map<std::string, std::string> subConfigMap;
    
    

public:
    
    int windowMaxSize = 100;
    int key = 0;
    int targetPartition = 0;  // The partition membership of a new vector, if the new vector is to be added to the window.
    std::vector<int> windowKeys;
    std::vector<int> partitionLabels;
    std::vector<int> nnIndices;  // A container to store the index of each point's nearest neighbor within the window.
    std::vector<double> nnDists;  // A container to store the nearest neighbor distance of each point within the window.
    std::unordered_map<int, double> avgNNDistPartitions;
    std::unordered_map<int, int> numPointsPartn;  // A dictionary to store the number of points in each partition.
    std::map<int, int> maxKeys;  // A dictionary to store the maxKey of each partition.
    std::vector<double> distsFromCurrVec;
    int keyToBeDeleted;
    int labelToBeDeleted;
    int indexToBeDeleted;
    double nnDistToBeDeleted;
    
    static pipePacket<nodeType> pPack;

    slidingWindow();
    pipePacket<nodeType> runPipe(pipePacket<nodeType>);
    void outputData(pipePacket<nodeType>);
    bool configPipe(std::map<std::string, std::string>);
    void runSubPipeline(pipePacket<nodeType>);
    void writeComplexStats(pipePacket<nodeType> &);
    static bool nnBasedEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
    static void deleteNNstats();
    static void updateStats();
};
