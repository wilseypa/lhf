#pragma once

// Header file for kd-Tree Pipe class - see kd-Tree.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class kdTreePipe : public basePipe {
  private:
	double beta;
  public:
    kdTreePipe();
    void runPipe(pipePacket& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket&);
};

