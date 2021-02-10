#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <vector>
#include <string>
#include "basePipe.hpp"
#include "simplexBase.hpp"
#include "utils.hpp"

class incrementalPersistence : public basePipe {
	private:
		int shift = 0;
		double maxEpsilon;
		std::map<std::string, std::string> configMap;
		std::string mode;

		struct cmpBySecond{ //Sort nodes by weight, then by lexicographic order
			bool operator()(simplexNode_P a, simplexNode_P b) const{
				if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
					return a->hash < b->hash;
				} else{
					return a->weight > b->weight;
				}
			}

			bool operator()(simplexNode* a, simplexNode* b) const{
				if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
					return a->hash < b->hash;
				} else{
					return a->weight > b->weight;
				}
			}
		};

	public:
		int dim;
	    incrementalPersistence();
	    void runPipe(pipePacket &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
		void outputData(pipePacket&);
};

