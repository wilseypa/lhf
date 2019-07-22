#pragma once

// Header file for basePipe class - see basePipe.cpp for descriptions
#include <map>
#include "pipePacket.hpp"
#include "utils.hpp"

class basePipe {
  private:
	utils ut;
  public:
	std::string pipeType;
	int debug;
    basePipe(){};
    basePipe(std::map<std::string, std::string> argDict);
    basePipe* newPipe(const std::string&);
    pipePacket runPipeWrapper(pipePacket);
	virtual void outputData(pipePacket);
    virtual pipePacket runPipe(pipePacket);
    virtual bool configPipe(std::map<std::string, std::string>);
};

