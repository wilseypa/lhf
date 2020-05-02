#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <vector>
#include "basePipe.hpp"
#include "simplexBase.hpp"

class unionFind{
	private:
		std::vector<int> rank, parent;
	public:
		unionFind(int n);
		int find(int i);
		bool join(int x, int y);
};

struct cmpBySecond{ //Sort nodes by weight, then by lexicographic order
	bool operator()(simplexBase::treeNode* a, simplexBase::treeNode* b) const{
		if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
			auto itA = a->simplex.rbegin(), itB = b->simplex.rbegin();
			while(itA != a->simplex.rend()){
				if(*itA != *itB) return *itA < *itB;
				++itA; ++itB;
			}
			return false;
		} else{
			return a->weight > b->weight;
		}
	}
};

class fastPersistence : public basePipe {
  private:
	utils ut;
	int shift = 0;
	double maxEpsilon;
  public:
	int dim;
    fastPersistence();
    pipePacket runPipe(pipePacket inData);
    bool configPipe(std::map<std::string, std::string> configMap);
	void outputData(pipePacket);
};

