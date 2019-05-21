#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <queue>
#include "basePipe.hpp"
#include "utils.hpp"
#include "indSimplexTree.hpp"


class optPersistencePairs : public basePipe {
  private:
  
	utils ut;
	double maxEpsilon;
	std::string twist;
  public:
	struct tArrayEntry_t{
		bool marked = false;
		int ti = -1;
		double birth;
		double death = -1;
		std::set<unsigned> simplex;
	};
	bool alterPipe = false;
	
	std::vector<tArrayEntry_t> tArray;
	
	int dim;
    optPersistencePairs();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	std::vector<std::vector<unsigned>> nSimplices(double, unsigned, std::vector<std::pair<double,std::vector<unsigned>>>);
	int checkFace(std::vector<unsigned> face, std::vector<unsigned>);
	int checkFace(std::set<unsigned> face, std::set<unsigned>);
	std::pair<std::set<unsigned>,std::set<unsigned>> getRankNull(std::vector<std::vector<unsigned>>);
	std::vector<std::vector<unsigned>> createBoundaryMatrix(std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges, int d, std::set<unsigned> pivots);
	std::vector<std::set<unsigned>> createBoundarySets(std::vector<std::vector<indSimplexTree::graphEntry>>, int, std::set<unsigned>, pipePacket);
	void outputData(pipePacket);
};

