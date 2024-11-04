#pragma once
#include <map>
#include <set>
#include <vector>
#include "helixPipe.hpp"

// base constructor
template <typename nodeType>
class helixDistPipe : public helixPipe<nodeType>
{
public:
  helixDistPipe();
  ~helixDistPipe();
  void runPipe(pipePacket<nodeType> &inData);
  bool configPipe(std::map<std::string, std::string> &configMap);
  void outputData(pipePacket<nodeType> &);
};
