#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "utils.hpp"


class bettiPipe : public basePipe {
  private:
	utils ut;
	float maxEpsilon;
  public:
	int dim;
    bettiPipe();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	std::vector<std::vector<unsigned>> nSimplices(double, unsigned, std::vector<std::pair<double,std::vector<unsigned>>>);
	int checkFace(std::vector<unsigned> face, std::vector<unsigned>);
	std::pair<int,int> reduceBoundaryMatrix(std::vector<std::vector<unsigned>>);
	std::pair<int,int> getRank(std::vector<std::vector<unsigned>>, std::vector<std::vector<unsigned>>);
};

