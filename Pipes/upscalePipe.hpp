#pragma once

// Header file for upscalePipe class - see upscalePipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class upscalePipe : public basePipe {
  private:
  public:
	int dim;
    upscalePipe();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
};

