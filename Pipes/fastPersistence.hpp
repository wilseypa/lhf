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

class binomialTable{
	private:
		std::vector<std::vector<long long>> v; 
	public:
		binomialTable(unsigned n, unsigned k);
		long long binom(unsigned n, unsigned k);
};

class fastPersistence : public basePipe {
  private:
	utils ut;
	int shift = 0;
	double maxEpsilon;
	long long ripsIndex(std::set<unsigned>& simplex, binomialTable& bin);
	unsigned maxVertex(long long ripsIndex, unsigned high, unsigned low, unsigned k, binomialTable &bin);
	std::vector<unsigned> getVertices(long long ripsIndex, int dim, unsigned n, binomialTable &bin);
  public:
	int dim;
    fastPersistence();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	void outputData(pipePacket);
};

