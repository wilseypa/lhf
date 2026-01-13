#pragma once

// Header file for betaPOlytopesPipe class - see betaPolytopes.cpp for descriptions
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <fstream>
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
	struct VectorHash;
	struct Simplex;
	struct Face;
	struct Strand;
	
	struct VectorHash {
    std::size_t operator()(const std::vector<unsigned>& v) const noexcept {
        std::size_t h = 0;
        for (unsigned x : v) {
            h ^= std::hash<unsigned>{}(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

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
		std::vector<std::shared_ptr<Simplex>> simplicies;

		Strand() = default;
	};

public:	
	betaPolytopes();
	void runPipe(pipePacket<nodeType> &inData);
	bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType> &);

	void flood_fill(Strand& strand, std::shared_ptr<Simplex>& simplex, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> &facelist);
	//helper functions for visualization:
	void collect_ids(const std::vector<Strand>& strands, std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids, std::unordered_map<const Simplex*, int>& simplex_ids);
	void write_faces_csv(const std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids);
	void write_simplices_csv(const std::vector<Strand>& strands, std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids, const std::unordered_map<const Simplex*, int>& simplex_ids);
	void export_strands_to_csv(const std::vector<Strand>& strands);


};
