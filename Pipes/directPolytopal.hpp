#pragma once

#include "basePipe.hpp"

// basePipe constructor
template <typename nodeType>
class directPolytopal : public basePipe<nodeType>
{
private:
public:
  directPolytopal();
  void runPipe(pipePacket<nodeType> &inData);
  bool configPipe(std::map<std::string, std::string> &configMap);
  void outputData(pipePacket<nodeType> &);
};
