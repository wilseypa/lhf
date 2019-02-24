#pragma once

// Header file for neighGraphPipe class - see neighGraphPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class neighGraphPipe : public basePipe {
  private:
	double epsilon;
  public:
    neighGraphPipe();
    pipePacket runPipe(pipePacket);
	void outputData(pipePacket);
    bool configPipe(std::map<std::string, std::string>);
};

