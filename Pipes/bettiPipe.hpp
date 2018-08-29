#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class bettiPipe : public basePipe {
  private:
  public:
	int dim;
    bettiPipe();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
};

