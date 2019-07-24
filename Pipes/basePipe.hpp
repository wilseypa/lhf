#pragma once

// Header file for basePipe class - see basePipe.cpp for descriptions
#include <map>
#include "pipePacket.hpp"
#include "utils.hpp"

class basePipe {
  private:
  public:
	utils ut;
	std::string pipeType;
	int debug;
	std::string outputFile;
    basePipe(){};
    basePipe* newPipe(const std::string&);
    pipePacket runPipeWrapper(pipePacket);
	virtual void outputData(pipePacket);
    virtual pipePacket runPipe(pipePacket);
    virtual bool configPipe(std::map<std::string, std::string>);
};

