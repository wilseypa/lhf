#pragma once

// Header file for slidingWindow class - see slidingWindow.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class slidingWindow : public basePipe {
  private:
	double epsilon;
	int dim;
	std::string inputFile;
	void runSubPipeline();
	std::map<std::string, std::string> subConfigMap;
  public:
    slidingWindow();
    pipePacket runPipe(pipePacket);
	void outputData(pipePacket);
    bool configPipe(std::map<std::string, std::string>);
    void runSubPipeline(pipePacket);
};

