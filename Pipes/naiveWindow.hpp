#pragma once

// Header file for naiveWindow class - see slidingWindow.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "utils.hpp"

template <class T>
class naiveWindow : public basePipe<T> {
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
    void runPipe(pipePacket<T>&);
    void outputData(pipePacket<T>&);
    bool configPipe(std::map<std::string, std::string>&);
    void runSubPipeline(pipePacket<T>);
    void writeComplexStats(pipePacket<T> &);
    static bool sampleStreamEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
};

template class naiveWindow<simplexNode>;
template class naiveWindow<alphaNode>;
