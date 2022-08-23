#pragma once

// Header file for naiveWindow class - see slidingWindow.cpp for descriptions
#include <map>
#include <unordered_map>
#include "basePipe.hpp"
#include "utils.hpp"
#define WINDOWNUMS 4 // constant number according to window number; 
#define INTERVAL 10
#define WINDOWSIZE 10
#define NEWDEBUGGER 0  //1 is cout or 0 is file output
#define PRCDEBUGGER true
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
        std::unordered_map<int, int> maxKeys;  //unordered version for maxkeys
        //std::map<int, int> maxKeys;  // A dictionary to store the maxKey of each partition.
        std::vector<double> distsFromCurrVec;
        int keyToBeDeleted;
        int labelToBeDeleted;
        int indexToBeDeleted;
        double nnDistToBeDeleted;
    };

    static EvalParams* defaultVals;
    static EvalParams* sumDefaultVals0; 
    static EvalParams* sumDefaultValsT[WINDOWNUMS+1]; 
    static EvalParams* totalDefaultVals;
    
    static pipePacket* pPack;
    static pipePacket* sumpPackT[WINDOWNUMS+1];
    static pipePacket* totalpPack;
    //static std::vector<std::vector<double>> sumoriginalData[WINDOWNUMS+1];
    
    static std::vector<std::vector<std::vector<double>>> sumoriginalData;
    
    slidingWindow();
    pipePacket runPipe(pipePacket);
    void outputData(pipePacket, int);
    void newoutputData(pipePacket, int);
    bool configPipe(std::map<std::string, std::string>);
    void runSubPipeline(pipePacket, int);
    void writeComplexStats(pipePacket &);
    static bool nnBasedEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
    static void deleteNNstats();
    static void updateStats();
    //===========================
    void test();
    void csyoutDebug(std::string, double, bool, bool );
    void csyoutDebug(std::string, int , bool, bool);
    void csyoutDebug(std::string, std::vector<double> , bool);
    void csyoutDebug(std::string, std::vector<int> , bool);
    void csyoutDebug(std::string, std::vector<int>, std::vector<double>, bool);
    void csyoutDebug(std::string, std::string, bool);
    void csyoutDebug(std::string, std::unordered_map<int, int>, bool);
    void csyoutDebug(std::string, std::unordered_map<int, double>, bool);
    void csyoutDebug(std::string,  std::vector<std::vector<double>> , bool);
    void csyoutDebugComplex(std::string, simplexBase*, bool);
    void computeTotal();
};

