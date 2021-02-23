#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <vector>
#include "basePipe.hpp"
#include "simplexBase.hpp"
#include "simplexArrayList.hpp"
#include "utils.hpp"

class involutedPersistence : public basePipe {
	private:
		std::vector<simplexNode*> edges;
		unsigned dim;

		struct sortLexicographic{ //Sort nodes by weight, then by lexicographic order
			bool operator()(simplexNode* a, simplexNode* b) const{
				if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
					auto itA = a->simplex.begin(), itB = b->simplex.begin();
					while(itA != a->simplex.end()){
						if(*itA != *itB) return *itA < *itB;
						++itA; ++itB;
					}
					return false;
				} else{
					return a->weight < b->weight;
				}
			}
		};

	public:
	    involutedPersistence();
	    void setupSimplices(std::vector<simplexNode*>, unsigned);
	    void runPipe(pipePacket &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
		void outputData(pipePacket&);
};

