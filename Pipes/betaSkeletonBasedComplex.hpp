#pragma once

// Header file for betaSkeletonBasedComplexPipe class - see betaSkeletonBasedComplexPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "../Preprocessing/kdTree.hpp"

class betaSkeletonBasedComplexPipe : public basePipe {
  private:
	double beta;
	double enclosingRadius;
	int dim;
	double epsilon;
  public:
    betaSkeletonBasedComplexPipe();
    void runPipe(pipePacket& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket&);
};

