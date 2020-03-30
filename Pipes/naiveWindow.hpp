#pragma once

// Header file for naiveWindow class - see slidingWindow.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "utils.hpp"

class naiveWindow : public basePipe {
  private:
	double epsilon;
	int dim;
	int repCounter = 0;
	std::string inputFile;
	void runSubPipeline();
	std::map<std::string, std::string> subConfigMap;
	void runComplexInitializer(pipePacket &);
  public:
    naiveWindow();
    pipePacket runPipe(pipePacket);
	void outputData(pipePacket);
    bool configPipe(std::map<std::string, std::string>);
    void runSubPipeline(pipePacket);
	void writeComplexStats(pipePacket &);
	static bool sampleStreamEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
};

