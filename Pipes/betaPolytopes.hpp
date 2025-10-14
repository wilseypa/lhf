#pragma once

// Header file for betaPOlytopesPipe class - see betaPolytopes.cpp for descriptions
#include <map>
#include <unordered_map>
#include <memory>
#include "basePipe.hpp"
#include "kdTree.hpp"

template <typename nodeType>
class betaPolytopes : public basePipe<nodeType>
{
private:
	double enclosingRadius;
	int dim;
	double epsilon;
	struct FacePtrHash;
	struct Simplex;
	struct Face;
	struct Strand;
	

	struct FacePtrHash {
    size_t operator()(const std::shared_ptr<Face>& f) const noexcept {
        std::hash<unsigned> hasher;
        size_t seed = 0;
        for (auto v : f->verticies) {
            seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    	}
	};

	struct FacePtrEq {
		bool operator()(const std::shared_ptr<Face>& a,
						const std::shared_ptr<Face>& b) const noexcept {
			return a->verticies == b->verticies;
		}
	};

	struct Simplex {
		std::vector<std::shared_ptr<Face>> faces;
		bool visited;

		Simplex() = default;
		Simplex(int d) : faces(), visited(false) {}
	};

	struct Face {
		std::vector<unsigned> verticies;
		std::vector<std::shared_ptr<Simplex>> adjacent_simplicies;

		bool operator==(const Face& other) const noexcept {
        return verticies == other.verticies;
		}

		Face() = default;
	};
	
	struct Strand {
		std::vector<Simplex> simplicies;

		Strand() = default;
	};

public:	
	betaPolytopes();
	void runPipe(pipePacket<nodeType> &inData);
	bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType> &);

	void flood_fill(Strand& strand, Simplex& simplex, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> &facelist);
};
