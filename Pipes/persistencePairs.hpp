#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <queue>
#include "basePipe.hpp"
#include "utils.hpp"


class persistencePairs : public basePipe {
  private:
	utils ut;
	float maxEpsilon;
	std::string twist;
  public:
	int dim;
    persistencePairs();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	std::vector<std::vector<unsigned>> nSimplices(double, unsigned, std::vector<std::pair<double,std::vector<unsigned>>>);
	int checkFace(std::vector<unsigned> face, std::vector<unsigned>);
	int checkFace(std::set<unsigned> face, std::set<unsigned>);
	void computeIntervals();
	void removePivotRows();
	void outputData(pipePacket);
};

