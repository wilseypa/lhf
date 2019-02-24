#pragma once

// Header file for ripsPipe class - see ripsPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class ripsPipe : public basePipe {
  private:
  public:
	int dim;
    ripsPipe();
    pipePacket runPipe(pipePacket);
    bool configPipe(std::map<std::string, std::string>);
	void outputData(pipePacket);
};

