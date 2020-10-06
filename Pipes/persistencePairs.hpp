#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <queue>
#include "basePipe.hpp"
#include "utils.hpp"


class persistencePairs : public basePipe {
  private:
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
    persistencePairs();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	std::vector<std::vector<unsigned>> nSimplices(double, unsigned, std::vector<std::pair<double,std::vector<unsigned>>>);
	int checkFace(std::vector<unsigned> face, std::vector<unsigned>);
	int checkFace(std::set<unsigned> face, std::set<unsigned>);
	std::set<unsigned> getRankNull(std::vector<std::vector<unsigned>>, unsigned);
	std::set<unsigned> createBoundaryMatrix(std::vector<std::vector<std::pair<std::set<unsigned>,double>>>);
	void outputData(pipePacket);
};

