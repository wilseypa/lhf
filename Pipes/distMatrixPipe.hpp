#pragma once

// Header file for distMatrixPipe class - see distMatrixPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class distMatrixPipe : public basePipe {
  private:
  public:
    distMatrixPipe();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	void outputData(pipePacket);
};

