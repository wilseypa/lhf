#pragma once

// Header file for streamVR class - see streamVR.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class streamVR : public basePipe {
  private:
	double epsilon;
	int dim;
  public:
    streamVR();
    pipePacket runPipe(pipePacket);
	void outputData(pipePacket);
    bool configPipe(std::map<std::string, std::string>);
};

