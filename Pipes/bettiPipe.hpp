#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <queue>
#include "basePipe.hpp"
#include "utils.hpp"


class bettiPipe : public basePipe {
  private:
	utils ut;
	float maxEpsilon;
	std::string twist;
  public:
	int dim;
    bettiPipe();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	std::vector<std::vector<unsigned>> nSimplices(double, unsigned, std::vector<std::pair<double,std::vector<unsigned>>>);
	int checkFace(std::vector<unsigned> face, std::vector<unsigned>);
	int checkFace(std::set<unsigned> face, std::set<unsigned>);
	std::pair<int,int> reduceBoundaryMatrix(std::vector<std::vector<unsigned>>);
	std::pair<std::queue<unsigned>, std::pair<int,int>> reduceBoundaryMatrixRev(std::vector<std::vector<unsigned>>);
	std::pair<std::queue<unsigned>, std::pair<int,int>> getRank(std::vector<std::vector<unsigned>>, std::vector<std::vector<unsigned>>, std::queue<unsigned>);
	std::pair<std::queue<unsigned>, std::pair<int,int>> removeSimplices(std::vector<std::pair<std::set<unsigned>,double>> &, std::vector<std::vector<unsigned>> &, double, std::queue<unsigned>&, std::pair<int, int>&);
	std::vector<std::vector<std::vector<unsigned>>> createBoundaryMatrix(std::vector<std::vector<std::pair<std::set<unsigned>,double>>>);
	void outputData(pipePacket);
};

