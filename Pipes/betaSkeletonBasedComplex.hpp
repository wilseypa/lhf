#pragma once

// Header file for betaSkeletonBasedComplexPipe class - see betaSkeletonBasedComplexPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "../Preprocessing/kdTree.hpp"

template <typename nodeType>
class betaSkeletonBasedComplexPipe : public basePipe<nodeType> {
  private:
	double beta;
	double enclosingRadius;
	int dim;
	double epsilon;
  public:
    betaSkeletonBasedComplexPipe();
    void runPipe(pipePacket<nodeType>& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType>&);
};

template class betaSkeletonBasedComplexPipe<simplexNode>;
template class betaSkeletonBasedComplexPipe<alphaNode>;
