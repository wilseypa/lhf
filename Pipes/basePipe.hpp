#pragma once

// Header file for basePipe class - see basePipe.cpp for descriptions
#include <map>
#include "pipePacket.hpp"

class basePipe {
  private:
  public:
	std::string pipeType;
	int debug;
    basePipe();
    basePipe* newPipe(const std::string&);
    pipePacket runPipeWrapper(pipePacket inData);
    virtual pipePacket runPipe(pipePacket inData);
    virtual bool configPipe(std::map<std::string, std::string> configMap);
};

