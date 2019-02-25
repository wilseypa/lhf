#pragma once

// Header file for boundaryPipe class - see boundaryPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "utils.hpp"

class boundaryPipe : public basePipe {
  private:
	std::vector<int> ranks;
	std::vector<int> nRanks;
	utils ut;
	float maxEpsilon;
  public:
	int dim;
    boundaryPipe();
    pipePacket runPipe(pipePacket);
    bool configPipe(std::map<std::string, std::string>);
	std::vector<std::vector<unsigned>> nSimplices(double, unsigned, std::vector<std::pair<double,std::vector<unsigned>>>);
	std::vector<std::vector<unsigned>> extractBoundaries(std::vector<std::pair<double, std::vector<unsigned>>>, std::vector<std::vector<unsigned>>, int);
	int checkFace(std::vector<unsigned> face, std::vector<unsigned>);
	void outputData(pipePacket);
	std::pair<std::vector<std::vector<unsigned>>,std::pair<int,int>> reduceBoundaryMatrix(std::vector<std::vector<unsigned>>);
	std::pair<std::vector<std::vector<unsigned>>,std::pair<int,int>> boundaryMatrix(std::vector<std::vector<unsigned>>, std::vector<std::vector<unsigned>>);
};

