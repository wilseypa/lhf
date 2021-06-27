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
	// For generation of combinations n choose r
	struct c_unique {
	unsigned current;
	c_unique() {current=0;}
	unsigned operator()() {return ++current;}
	} UniqueNumber;
  public:
    betaSkeletonBasedComplexPipe();
    void runPipe(pipePacket<nodeType>& inData);
    bool checkInsertDsimplex(std::vector<unsigned> dsimplex,pipePacket<nodeType> &inData,double beta,double averageDistance,kdTree tree);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType>&);
};
