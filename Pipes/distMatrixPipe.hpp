#pragma once

// Header file for distMatrixPipe class - see distMatrixPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

template<typename nodeType>
class distMatrixPipe : public basePipe<nodeType> {
  private:
	double enclosingRadius;
  public:
    distMatrixPipe();
    void runPipe(pipePacket<nodeType>& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType>&);
};


//Explicit Template Class Instantiation
template class distMatrixPipe<simplexNode>;
template class distMatrixPipe<alphaNode>;
