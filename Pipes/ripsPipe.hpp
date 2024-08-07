#pragma once

// Header file for ripsPipe class - see ripsPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

template <typename nodeType>
class ripsPipe : public basePipe<nodeType>
{
private:
public:
  std::string collapse;
  int dim;
  ripsPipe();
  void runPipe(pipePacket<nodeType> &);
  bool configPipe(std::map<std::string, std::string> &);
  void outputData(pipePacket<nodeType> &);
};
