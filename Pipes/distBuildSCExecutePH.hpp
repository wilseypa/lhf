#pragma once

// Header file for buldingValidSimplicialComplexes class - see buildComplexes.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class distBuildSCExecutePHPipe : public basePipe {
	private:
 	public:
		int dim;
	    distBuildSCExecutePHPipe();
		std::vector<set<simplexNode_P, cmpByWeight>> distBuildSCExecutePHPipe::buildValidSimplicialComplex(vector<set<unsigned>> dsimplexes,pipePacket &inData){
        pair<T, bool> getNthElement(set<T>& searchSet,unsigned index);
	    void runPipe(pipePacket &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
		void outputData(pipePacket&);
};
