#pragma once

// Header file for betaPOlytopesPipe class - see betaPolytopes.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "kdTree.hpp"

template <typename nodeType>
class betaPolytopes : public basePipe<nodeType>
{
private:
	double enclosingRadius;
	int dim;
	double epsilon;

	struct VectorHash {
		size_t operator()(const std::vector<unsigned>& v) const noexcept {
			std::hash<int> hasher;
			size_t seed = 0;
			for (int i : v) {
				seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};
	
	struct Simplex {
		std::vector<unsigned> simplex;
		std::vector<std::vector<unsigned>> faces;
		bool visited;

		Simplex(int d) : simplex(d+1), faces(d+1, std::vector<unsigned>(d)), visited(false) {}

	};

public:	
	betaPolytopes();
	void runPipe(pipePacket<nodeType> &inData);
	bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType> &);
};
