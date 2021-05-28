#pragma once

// Header file for betaSkeletonBasedComplexPipe class - see betaSkeletonBasedComplexPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "../Preprocessing/kdTree.hpp"

template <typename T>
class betaSkeletonBasedComplexPipe : public basePipe<T> {
  private:
	double beta;
	double enclosingRadius;
	int dim;
	double epsilon;
  public:
    betaSkeletonBasedComplexPipe();
    void runPipe(pipePacket<T>& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<T>&);
};

template class betaSkeletonBasedComplexPipe<simplexNode>;
template class betaSkeletonBasedComplexPipe<alphaNode>;
