#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <vector>
#include "basePipe.hpp"

class unionFind{
	private:
		std::vector<int> rank, parent;
	public:
		unionFind(int n);
		int find(int i);
		bool join(int x, int y);
};

class fastPersistence : public basePipe {
  private:
	utils ut;
	double maxEpsilon;
  public:		
	int dim;
    fastPersistence();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	void outputData(pipePacket);
};

