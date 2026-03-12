#pragma once

// Header file for betaPOlytopesPipe class - see betaPolytopes.cpp for descriptions
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <fstream>
#include <vector>
#include <queue>
#include <set>
#include <Eigen/Dense>
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/Qhullfacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
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
        for (auto v : f->vertices) {
            seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    	}
	};

	struct FacePtrEq {
		bool operator()(const std::shared_ptr<Face>& a,
						const std::shared_ptr<Face>& b) const noexcept {
			return a->vertices == b->vertices;
		}
	};

	struct Simplex {
		std::vector<std::shared_ptr<Face>> faces;
		bool visited;

		Simplex() = default;
		Simplex(int d) : faces(), visited(false) {}
	};

	struct Face {
		std::vector<unsigned> vertices;
		std::vector<std::shared_ptr<Simplex>> adjacent_simplices;

		bool operator==(const Face& other) const noexcept {
        return vertices == other.vertices;
		}

		Face() = default;
	};
	
	struct Strand {
		std::vector<std::shared_ptr<Simplex>> simplices;
    	std::vector<typename betaPolytopes<nodeType>::Chart> atlas;

		Strand() = default;
	};

	struct Polytope {
		std::vector<Eigen::VectorXd> vertices;
		std::vector<std::vector<int>> faces;
	};

	struct Chart {
		int intrinsic_dim; //intrinsic dim
		int d; //ambient dim
		size_t num_points; //num of point accumulated
		double distortion_threshold;
		Polytope polytope;

		int global_id = -1;
		int strand_id = -1;

		//running PCA state
		Eigen::VectorXd mean; //dx1
		Eigen::MatrixXd M2; //dxd

		//current basis
		Eigen::MatrixXd basis;

		//assigned simplices
		std::vector<std::shared_ptr<Simplex>> simplices;

		Chart(int intrinsic_dim, int ambient_dim, double threshold)
			: intrinsic_dim(intrinsic_dim),
			  d(ambient_dim),
			  num_points(0),
			  distortion_threshold(threshold),
			  mean(Eigen::VectorXd::Zero(ambient_dim)),
			  M2(Eigen::MatrixXd::Zero(ambient_dim, ambient_dim)),
			  basis(Eigen::MatrixXd::Zero(ambient_dim, intrinsic_dim))
		{}	

		bool tryAddSimplex(const std::shared_ptr<Simplex>& simplex, const std::vector<Eigen::VectorXd>& cloud_points);

	};

public:	
	betaPolytopes();
	void runPipe(pipePacket<nodeType> &inData);
	bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType> &);

	void flood_fill(Strand& strand, std::shared_ptr<Simplex>& simplex, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> &facelist);	
	std::vector<typename betaPolytopes<nodeType>::Chart> generateAtlasForStrand(std::vector<std::shared_ptr<Simplex>> mesh_structs, const std::vector<Eigen::VectorXd>& cloud_points, int intrinsic_dim, double distortion_threshold, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> &facelist);
	void assignChartIds(std::vector<Strand>& strands);
	std::vector<int> collectChartVertexIndices(const Chart& chart);
	bool canMerge(Chart& A, Chart& B, const std::vector<Eigen::VectorXd>& cloud_points);
	Polytope computeConvexHull(const std::vector<Eigen::VectorXd>& cloud_points, const std::vector<int>& vertex_indices);
	void computeHullForChart(Chart& chart, const std::vector<Eigen::VectorXd>& cloud_points);

	//helper functions for visualization:
	void exportAtlasStructure(const std::vector<Chart>& atlas);
	void collect_ids(const std::vector<Strand>& strands, std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids, std::unordered_map<const Simplex*, int>& simplex_ids);
	void write_faces_csv(const std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids);
	void write_simplices_csv(const std::vector<Strand>& strands, std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids, const std::unordered_map<const Simplex*, int>& simplex_ids);
	void export_strands_to_csv(const std::vector<Strand>& strands);


};
