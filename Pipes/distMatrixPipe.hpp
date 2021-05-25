#pragma once

// Header file for distMatrixPipe class - see distMatrixPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

template<typename T>
class distMatrixPipe : public basePipe<T> {
  private:
	double enclosingRadius;
  public:
    distMatrixPipe();
    void runPipe(pipePacket<T>& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<T>&);
};


//Explicit Template Class Instantiation
template class distMatrixPipe<simplexNode>;
template class distMatrixPipe<alphaNode>;
