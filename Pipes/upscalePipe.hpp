#pragma once

// Header file for upscalePipe class - see upscalePipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

template <typename T>
class upscalePipe : public basePipe<T> {
  private:
  public:
	std::map<std::string, std::string> subConfigMap;
	int dim;
	double scalarV;
    upscalePipe();
	void runSubPipeline(pipePacket<T>&);
	//void outputData(pipePacket<T>&);
    void runPipe(pipePacket<T>&);
    bool configPipe(std::map<std::string, std::string>&);    
};

template class upscalePipe<simplexNode>;
template class upscalePipe<alphaNode>;

