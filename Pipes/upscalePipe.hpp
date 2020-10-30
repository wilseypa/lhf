#pragma once

// Header file for upscalePipe class - see upscalePipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class upscalePipe : public basePipe {
  private:
  public:
	std::map<std::string, std::string> subConfigMap;
	int dim;
	double scalarV;
    upscalePipe();
	void runSubPipeline(pipePacket&);
	//void outputData(pipePacket&);
    void runPipe(pipePacket&);
    bool configPipe(std::map<std::string, std::string>&);    
};

