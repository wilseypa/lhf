#pragma once

// Header file for upscalePipe class - see upscalePipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

template <typename nodeType>
class upscalePipe : public basePipe<nodeType> {
  private:
  public:
	std::map<std::string, std::string> subConfigMap;
	int dim;
    int maxSize;
	double scalarV;
    upscalePipe();
	void runSubPipeline(pipePacket<nodeType>&);
	//void outputData(pipePacket<nodeType>&);
    void runPipe(pipePacket<nodeType>&);
    bool configPipe(std::map<std::string, std::string>&);    
};
