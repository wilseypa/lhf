#pragma once

// Header file for betaSubSkeletonComplexPipe class - see betaSubSkeletonComplexPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "kdTree.hpp"

template <typename nodeType>
class betaSubSkeletonComplex : public basePipe<nodeType>
{
private:
	double beta;
	std::string betaMode;
	double enclosingRadius;
	int dim;
	std::string betaMesh;
	double epsilon;
	// For generation of combinations n choose r
	struct c_unique
	{
		unsigned current;
		c_unique() { current = 0; }
		unsigned operator()() { return ++current; }
	} UniqueNumber;

public:
	betaSubSkeletonComplex();
	void runPipe(pipePacket<nodeType> &inData);
	bool checkInsertSubDsimplex(std::vector<unsigned> dsimplex, pipePacket<nodeType> &inData, double beta, double averageDistance, kdTree tree);
	bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType> &);
};
