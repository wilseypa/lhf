#pragma once

// Header file for naiveWindow class - see slidingWindow.cpp for descriptions
#include <map>
#include <unordered_map>
#include "basePipe.hpp"
#include "utils.hpp"

class slidingWindow : public basePipe {
private:
    double epsilon;
    int dim;
    int repCounter = 0;
    std::string inputFile;
    void runSubPipeline();
    std::map<std::string, std::string> subConfigMap;

public:
    struct EvalParams
    {
        int windowMaxSize;
        int key;
        int targetPartition;  // The partition membership of a new vector, if the new vector is to be added to the window.
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
    };

    static EvalParams* defaultVals;
    static pipePacket* pPack;

    slidingWindow();
    pipePacket runPipe(pipePacket);
    void outputData(pipePacket);
    bool configPipe(std::map<std::string, std::string>);
    void runSubPipeline(pipePacket);
    void writeComplexStats(pipePacket &);
    static bool nnBasedEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
    static void deleteNNstats();
    static void updateStats();
};

