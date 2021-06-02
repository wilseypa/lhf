#pragma once

// Header file for neighGraphPipe class - see neighGraphPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

template<typename nodeType>
class neighGraphPipe : public basePipe<nodeType> {
  private:
  	double epsilon;
	  int dim;
  public:
    neighGraphPipe();
    void runPipe(pipePacket<nodeType>&);
	  void outputData(pipePacket<nodeType>&);
    bool configPipe(std::map<std::string, std::string>&);
};
