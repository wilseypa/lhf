#pragma once

// Header file for ripsPipe class - see ripsPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"


template<typename T>
class ripsPipe : public basePipe<T> {
  private:
  public:
    std::string collapse;
    int dim;
    ripsPipe();
    void runPipe(pipePacket<T>&);
    bool configPipe(std::map<std::string, std::string>&);
    void outputData(pipePacket<T>&);
};

template class ripsPipe<simplexNode>;
template class ripsPipe<alphaNode>;
