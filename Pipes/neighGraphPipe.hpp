#pragma once

// Header file for neighGraphPipe class - see neighGraphPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

template<typename T>
class neighGraphPipe : public basePipe<T> {
  private:
	double epsilon;
	int dim;
  public:
    neighGraphPipe();
    void runPipe(pipePacket<T>&);
	void outputData(pipePacket<T>&);
    bool configPipe(std::map<std::string, std::string>&);
};


//Explicit Template Class Instantiation
template class neighGraphPipe<simplexNode>;
template class neighGraphPipe<alphaNode>;
