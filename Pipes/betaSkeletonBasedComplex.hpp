#pragma once

// Header file for betaSkeletonBasedComplexPipe class - see betaSkeletonBasedComplexPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "../Preprocessing/kdTree.hpp"

template <typename nodeType>
class betaSkeletonBasedComplex : public basePipe<nodeType> {
  private:
	double beta;
	std::string betaMode;
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
    betaSkeletonBasedComplex();
    void runPipe(pipePacket<nodeType>& inData);
    bool checkInsertDsimplex(std::vector<unsigned> dsimplex,pipePacket<nodeType> &inData,double beta,double averageDistance,kdTree tree);
    unsigned selectCenter(std::vector<double> hpcofffaces, std::vector<std::vector<double>> betaCenters,std::vector<double> otherPoint);
    bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType>&);
};
