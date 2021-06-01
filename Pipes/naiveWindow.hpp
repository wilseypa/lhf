#pragma once

// Header file for naiveWindow class - see slidingWindow.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "utils.hpp"

template <class nodeType>
class naiveWindow : public basePipe<nodeType> {
  private:
	double epsilon;
	int dim;
	int repCounter = 0;
	std::string inputFile;
	void runSubPipeline();
	std::map<std::string, std::string> subConfigMap;
	//void runComplexInitializer(pipePacket &);
  public:
    std::vector<std::vector<double>> distMatrix;
    naiveWindow();
    void runPipe(pipePacket<nodeType>&);
    void outputData(pipePacket<nodeType>&);
    bool configPipe(std::map<std::string, std::string>&);
    void runSubPipeline(pipePacket<nodeType>);
    void writeComplexStats(pipePacket<nodeType> &);
    static bool sampleStreamEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
};

template class naiveWindow<simplexNode>;
template class naiveWindow<alphaNode>;
