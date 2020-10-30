#pragma once

// Header file for distMatrixPipe class - see distMatrixPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class distMatrixPipe : public basePipe {
  private:
	double enclosingRadius;
  public:
    distMatrixPipe();
    void runPipe(pipePacket& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket&);
};

