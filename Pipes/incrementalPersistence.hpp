#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <vector>
#include <string>
#include "basePipe.hpp"
#include "simplexBase.hpp"
#include "utils.hpp"

template <typename nodeType>
class incrementalPersistence : public basePipe<nodeType> {
	private:
		int shift = 0;
		unsigned nPts = 0;
		double maxEpsilon;
		std::string inv = "false";

		bool saveVertices = false; //Should we save the vertices of the simplices, or just their hashes

		struct sortReverseLexicographic{ //Sort nodes by weight, then by reverse lexicographic order
			template <class simplexNodePointer>
			bool operator()(simplexNodePointer a, simplexNodePointer b) const{
				if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order
					return a->hash < b->hash;
				} else{
					return a->weight > b->weight;
				}
			}
		};

		struct sortLexicographic{ //Sort nodes by weight, then by lexicographic order
			template <class simplexNodePointer>
			bool operator()(simplexNodePointer a, simplexNodePointer b) const{
				if(a->weight == b->weight){ //If the simplices have the same weight, sort by lexicographic order
					return a->hash < b->hash;
				} else{
					return a->weight < b->weight;
				}
			}
		};
		
		
		typedef std::shared_ptr<nodeType> templateNode_P;

	public:
		int dim;
	    incrementalPersistence();
	    void runPipe(pipePacket<nodeType> &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
		void outputData(pipePacket<nodeType>&);

		template <typename simplexNodePointer, typename comp>
		std::vector<simplexNodePointer> incrementalByDimension(pipePacket<nodeType>&, std::vector<simplexNodePointer>&, std::vector<simplexNodePointer> pivots, unsigned, comp, std::string, bool);
};
