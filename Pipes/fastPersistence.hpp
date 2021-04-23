#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <vector>
#include "basePipe.hpp"
#include "simplexBase.hpp"
#include "utils.hpp"

class fastPersistence : public basePipe {
	private:
		int shift = 0;
		double maxEpsilon;
		std::string inv = "false";

		struct sortReverseLexicographic{ //Sort nodes by weight, then by lexicographic order
			bool operator()(simplexNode_P a, simplexNode_P b) const{
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
		
		struct sortLexicographic{ //Sort nodes by weight, then by lexicographic order
			bool operator()(simplexNode_P a, simplexNode_P b) const{
				if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
					auto itA = a->simplex.rbegin(), itB = b->simplex.rbegin();
					while(itA != a->simplex.rend()){
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
		int dim;
	    fastPersistence();
	    void runPipe(pipePacket &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
		void outputData(pipePacket&);

		template <class simplexNodePointer, class comp>
		std::vector<simplexNodePointer> persistenceByDimension(pipePacket&, std::vector<simplexNodePointer>, std::vector<simplexNodePointer> pivots, unsigned, comp, std::string, bool);
};
