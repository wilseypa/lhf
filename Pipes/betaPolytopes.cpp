	/*
 * betaSubSkeletonComplexPipe hpp + cpp extend the basePipe class for calculating the
 * beta Skeleton Based Complex generation for data input
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <set>
#include <algorithm>
#include <unordered_map>
#include "betaPolytopes.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include "readInput.hpp"

using orgQhull::Qhull;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexSet;

// basePipe constructor
template <typename nodeType>
betaPolytopes<nodeType>::betaPolytopes()
{
	this->pipeType = "betaPolytopes";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void betaPolytopes<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	std::vector<std::vector<unsigned>> dsimplexmesh = inData.dsimplexmesh;
	std::vector<std::shared_ptr<Simplex>> mesh_structs; //copy in here to have struct features in each simplex
	std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> facelist; //stores <face, #of incident simplex>
	std::vector<Strand> strands;
	std::cout << "Simplices in beta mesh::" << dsimplexmesh.size()<<std::endl;

	std::ofstream out("betaMesh.csv");
	int counter = 0;


	for (auto& row : dsimplexmesh) {
		for (auto col : row)
			out << col << ',';
		out << '\n';
	}

	for(auto x:dsimplexmesh){
		counter++;
	}
	std::cout << counter << std::endl;

	std::ofstream outFile("../../python_tests/Polytopal_Development/vertices.csv");

// Check if the file opened successfully
	if (!outFile) {
		std::cerr << "Error opening vertices.csv." << std::endl;
		return;
	}

	// Write CSV header
	outFile << "vertex_id";
	if (!dsimplexmesh.empty()) {
		for (size_t d = 0; d < dsimplexmesh[0].size(); ++d) {
			outFile << ",c" << d;
		}
	}
	outFile << "\n";

	// Write vertex rows
	for (size_t vid = 0; vid < dsimplexmesh.size(); ++vid) {
		outFile << vid;
		for (const auto& coord : dsimplexmesh[vid]) {
			outFile << "," << coord;
		}
		outFile << "\n";
	}

	outFile.close();

	
	std::cout<<"We will generate Polytopes here"<<std::endl;
	
	//iterate through each simplex and construct each Simplex object as list of its faces
	for(const auto& x:dsimplexmesh) {
		auto current_simplex = std::make_shared<Simplex>(this->dim);
		//iterate through each face of each simplex to construct Face object as vertices and adjacent simplices
		for(int i = 0; i < x.size(); ++i) {
			auto current_face = std::make_shared<Face>();
			for(int j = 0; j < x.size(); ++j) {
				if (i != j){
					current_face -> vertices.push_back(x[j]);
				}
			}
			current_simplex -> faces.push_back(current_face);

			auto [it, inserted] = facelist.emplace(current_face, 1); //facelist as hashtable makes this check fast
			if(!inserted) {
				it -> first -> adjacent_simplices.push_back(current_simplex);
				it -> second++;
			}
			else {current_face -> adjacent_simplices.push_back(current_simplex);}
		}
		mesh_structs.push_back(current_simplex);
	}

	//convert input data to eigenvecs	
	std::vector<Eigen::VectorXd> cloud_points;
	cloud_points.reserve(inData.inputData.size());

	for (const auto& row : inData.inputData) {
		Eigen::VectorXd vec = Eigen::Map<const Eigen::VectorXd>(row.data(), row.size());
		cloud_points.push_back(vec);
	}

	std::vector<typename betaPolytopes<nodeType>::Chart> atlas = generateAtlasForStrand(mesh_structs, cloud_points, 2, 0.05, facelist);

	exportAtlasStructure(atlas);
	/*
	//We will still leverage this section to identify the structure of the manifolds in the mesh
	// strand enumeration code, reviewing necessity in algorithm
	for(auto simplex:mesh_structs) {
		if(simplex -> visited == false) {
			Strand strand; //whenever last recursion ends (last strand fully enumerated), find a new unvisited simplex to enumerate new strand
			flood_fill(strand, simplex, facelist);
			strands.push_back(strand);
		}
	}
	

	export_strands_to_csv(strands);

	*/

	//std::cout << "Generating atlases from strands..." << std::endl;
/*
	//convert input data to eigenvecs	
	std::vector<Eigen::VectorXd> cloud_points;
	cloud_points.reserve(inData.inputData.size());

	for (const auto& row : inData.inputData) {
		Eigen::VectorXd vec = Eigen::Map<const Eigen::VectorXd>(row.data(), row.size());
		cloud_points.push_back(vec);
	}

	//iterate over each strand and build atlas of non-overlapping charts
	//Will want to make ambient dim = ambient dim of structure.
	//can consider making distortion factor configurable from cmd line.
	
	for(auto& strand : strands) {
		generateAtlasForStrand(strand, cloud_points, 2, 0.1, facelist);
	}

	//assigning ids and exporting atlases to files
	assignChartIds(strands);

	// build simplex -> id map
	std::unordered_map<const Simplex*, int> simplex_ids;
	int next_simplex_id = 0;

	for (const auto& strand : strands) {
		for (const auto& simplex : strand.simplices) {
			if (!simplex_ids.count(simplex.get())) {
				simplex_ids[simplex.get()] = next_simplex_id++;
			}
		}
	}
	

	std::ofstream atlasFile("../../python_tests/Polytopal_Development/atlases.csv");
	atlasFile << "strand_id,chart_id,simplex_id\n";

	for (size_t strand_id = 0; strand_id < strands.size(); ++strand_id) {
		const auto& strand = strands[strand_id];

		for (const auto& chart : strand.atlas) {
			for (const auto& simplex : chart.simplices) {
				int sid = simplex_ids[simplex.get()];

				atlasFile << chart.strand_id << "," << chart.global_id << "," << sid << "\n";
			}
		}
	}

	atlasFile.close();

	std::cout << "Atlas CSV written successfully." << std::endl;

	//build mapping of which simplex each chart belongs to
	//O(n) step with respect to num simplices in mesh
	std::unordered_map<const Simplex*, int> simplexToChart;

	for (const auto& strand : strands) {
		for (const auto& chart : strand.atlas) {
			for (const auto& simplex : chart.simplices) {
				simplexToChart[simplex.get()] = chart.global_id;
			}
		}
	}

	//build chart adjacency graph (only across strands)
	std::unordered_map<int, std::unordered_set<int>> chartAdj;

	std::unordered_map<int, int> chartToStrand; //chart to strand mapping
	for (const auto& strand : strands) {
		for (const auto& chart : strand.atlas) {
			chartToStrand[chart.global_id] = chart.strand_id;
		}
	}

	//O(Num of simplices * faces per simplex) linear behavior for mesh
	for (const auto& strand : strands) {
		for (const auto& chart : strand.atlas) {
			int chart_id = chart.global_id;

			for (const auto& simplex : chart.simplices) {
				for (const auto& face : simplex -> faces) {
					for (const auto& neighbor_simplex : face -> adjacent_simplices) {
						const Simplex* neighbor_ptr = neighbor_simplex.get();
						if (neighbor_ptr == simplex.get()) continue; //skip self

						int neighbor_chart_id = simplexToChart[neighbor_ptr];

						//only consider different charts
						if (neighbor_chart_id == chart_id) continue;

						//only across strands
						if (chartToStrand[neighbor_chart_id] == chartToStrand[chart_id]) continue;

						//add adjacency
						chartAdj[chart_id].insert(neighbor_chart_id);
						chartAdj[neighbor_chart_id].insert(chart_id);
					}
				}
			}
		}
	}

	*/


//TESTS
/*
	std::cout << "Printing strands:\n";

	int strand_idx = 0;
	for (const auto& strand : strands) {
		std::cout << "Strand " << strand_idx++ << ":\n";

		int simplex_idx = 0;
		for (const auto& simplex : strand.simplices) {
			std::cout << "  Simplex " << simplex_idx++ << ":\n";

			int face_idx = 0;
			for (const auto& face : simplex->faces) {
				std::cout << "    Face " << face_idx++ << ": ";

				for (const auto& vert : face->vertices) {
					std::cout << vert << " ";
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
*/
	
	/* Outlie of the algorithm that I have in mind.
	1. Intialize every simplex in the mesh as unvisited
	2. Start from a random simplex which is unvisited
	    3. Initialize a strand as current simplex and strand boundary as simplex facets
	      4. For every facet f in the boundary.
				remove f from boundary
	            5. Find the count(n) of cofacets (cfs) of facet f other than those in the strand
				6.   if(n==1):
				         add that cofacet to the strand and its unexplored facets to boundary.
				     else if(n>1)
						 This is bifurcation, trifucrcation etc. stop. The stand will not grow in that directon
	
	After strands are identified it is time to flatten them and do convex decomposition.
	
	Outline of changes so far after strand generation
	Each strand is identified as an atlas of 1 or more non-overlapping charts
	Current approach does the following per strand:
		init a chart, add neighboring simplices to a queue
		for each simplex in queue, test distortion when adding to chart, reject if too large
		if we add to chart, add neighboring simplices to queue for potential chart enumeration (bfs essentially)
	Each strand is considered an atlas constructed of charts with the above ideas

	Idea right now is each chart is to be decomposed into a polytope.
	Idea of next steps are as follows:
		for each chart, collect all points
		use the PCA basis for the chart to project into intrinsic dim
		compute convex hull of points in that space
		build the polytope

	I will implement this soon + python scripts to visualize, we can then begin computing homologies and ironing out potential issues.
		
	The strands seem to be nearly 1:1 with simplices in the case of a 2d mesh. What's interesting is I noticed it grows
	quite a bit with higher dimensional meshes. When we utilize higher d meshes, we gain accuracy with cost 
	of higher complexity. The strand -> chart -> poly pipeline might be better served for higher d mesh cases
	wondering if we should configure such that we jump from mesh -> convex hull calculations in the 2d mesh case
	*/
	
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool betaPolytopes<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
{
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if (pipe != configMap.end())
	{
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if (pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();

	this->ut = utils(strDebug, this->outputFile);
	pipe = configMap.find("dimensions");
	if (pipe != configMap.end())
	{
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else
		return false;

	this->configured = true;
	this->ut.writeDebug("betaPolytopes Pipe ", "Configured with parameters { eps: " + configMap["epsilon"] + configMap["beta"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void betaPolytopes<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	// Output related to betaPolytopes
	return;
}

template <typename nodeType>
void betaPolytopes<nodeType>::flood_fill(Strand& strand, std::shared_ptr<Simplex>& simplex, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq>& facelist) {
	if(simplex -> visited) return; //avoid duplicate simplices
	strand.simplices.push_back(simplex);
	simplex -> visited = true; //mark as visited once pushed to a strand
	//for each face, check if not an intersection, then recursively call function to add adjacent simplex to strand
	for(const auto& face:simplex -> faces) {
		if(facelist.at(face)==2) {
			for(auto& current_simplex:face -> adjacent_simplices) {
				flood_fill(strand, current_simplex, facelist);
			}
		}
	}
	return;
}

template <typename nodeType>
bool betaPolytopes<nodeType>::Chart::tryAddSimplex(const std::shared_ptr<Simplex>& simplex, const std::vector<Eigen::VectorXd>& cloud_points) {
	//collect simplex points
	std::unordered_set<unsigned> vertex_ids;
	for (auto face : simplex->faces) {
		for (auto vid : face->vertices) {
			vertex_ids.insert(vid);
		}
	}
	std::vector<Eigen::VectorXd> points;
	for (auto vid : vertex_ids) {
		points.push_back(cloud_points[vid]);
	}
	
	//accept if chart is empty automatically
	if (num_points == 0) {
		simplices.push_back(simplex);
		for (const auto& x : points) {
			num_points++;
			Eigen::VectorXd delta = x - mean;
			mean += delta / num_points;
			Eigen::VectorXd delta2 = x - mean;
			M2 += delta * delta2.transpose();
		}
		return true;
	}

	//test distortion against current chart basis
	for (const auto& x : points) {
		Eigen::VectorXd x_c = x - mean;
		Eigen::VectorXd proj = basis * basis.transpose() * x_c;
		double dist = (x_c - proj).norm();
		std::cout << "distortion " << dist << "\n";

		if (dist > distortion_threshold) {
			return false; //reject simplex
		}
	}

	//distortion test passed, update chart
	Eigen::VectorXd new_mean = mean;
	size_t total_points = num_points + points.size();

	for (const auto& x : points) {
		new_mean += (x - new_mean) / total_points;
	}

	Eigen::MatrixXd M2_new = M2;
	for (const auto& x : points) {
		Eigen::VectorXd delta1 = x - mean;
		Eigen::VectorXd delta2 = x - new_mean;
		M2_new += delta1 * delta2.transpose();
	}

	//recompute PCA basis
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M2_new / (total_points - 1));
	Eigen::MatrixXd U = solver.eigenvectors().rightCols(intrinsic_dim);

	//accept simplex
	simplices.push_back(simplex);
	mean = new_mean;
	M2 = M2_new;
	basis = U;
	num_points = total_points;

	return true;
}

//currently we measure distortion on a per simplex basis. If we have accuracy issues, one area to look
//will be monitoring global distortion to avoid poorly conditioned charts for large strands
template <typename nodeType>
std::vector<typename betaPolytopes<nodeType>::Chart> betaPolytopes<nodeType>::generateAtlasForStrand(std::vector<std::shared_ptr<Simplex>> mesh_structs, const std::vector<Eigen::VectorXd>& cloud_points, int intrinsic_dim, double distortion_threshold, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> &facelist) {
	std::vector<typename betaPolytopes<nodeType>::Chart> atlas;
	
	//reset simplices visited
	for(auto& simplex : mesh_structs) {
		simplex->visited = false;
	}

	//iterate through until all simplices in strand are assigned
	for (auto& seed : mesh_structs) {
		if (seed -> visited) continue;

		//create new chart
		Chart chart(seed->faces[0]->vertices.size(), intrinsic_dim, distortion_threshold);

		//init chart
		if(!chart.tryAddSimplex(seed, cloud_points)) {
			//if this actually happens distortion threshold may be too small
			seed->visited = true;
			continue;
		}
		seed->visited = true;

		//bfs queue
		std::queue<std::shared_ptr<Simplex>> q;
		q.push(seed);

		while(!q.empty()) {
			auto current = q.front();
			q.pop();

			//get neighbors. enforces enumeration on manifolds, we can actually remove strand enumeration just using this check instead.
			for(const auto& face:current -> faces) {
				//if(facelist.at(face) != 2) continue;

				for(const auto& neighbor : face->adjacent_simplices) {

					if(neighbor == current) continue;

					if(neighbor->visited) continue;

					if(chart.tryAddSimplex(neighbor, cloud_points)) {
						neighbor->visited = true;
						q.push(neighbor);
					}
				}
			}
		}
		//store chart in atlas
		atlas.push_back(chart);
	}
	return atlas;
}

template <typename nodeType>
void betaPolytopes<nodeType>::assignChartIds(std::vector<Strand>& strands) {
	int next_id = 0;

	for (int s = 0; s < strands.size(); ++s) {
		Strand& strand = strands[s];

		for (Chart& chart : strand.atlas) {
			chart.global_id = next_id;
			chart.strand_id = s;
			next_id++;
		}
	}
}

template <typename nodeType>
std::vector<int> betaPolytopes<nodeType>::collectChartVertexIndices(const Chart& chart) {
	std::unordered_set<int> unique_ids;

	for (const std::shared_ptr<Simplex>& simplex : chart.simplices) {
		for (const std::shared_ptr<Face>& face : simplex->faces) {
			for (int vid : face->vertices) {
				unique_ids.insert(vid);
			}
		}
	}
	return std::vector<int>(unique_ids.begin(), unique_ids.end());
}

template <typename nodeType>
typename betaPolytopes<nodeType>::Polytope betaPolytopes<nodeType>::computeConvexHull(const std::vector<Eigen::VectorXd>& cloud_points, const std::vector<int>& vertex_indices) {
	Polytope poly;

	if (vertex_indices.empty())
		return poly;
	
	const int dim = cloud_points[0].size();
	const int num_points = vertex_indices.size();

	// flatten points into contigous buffer for Qhull
	std::vector<double> coords;
	coords.reserve(num_points * dim);

	for (int vid : vertex_indices) {
		for (int d = 0; d < dim; ++d) {
			coords.push_back(cloud_points[vid](d)); //ensure correctness here
		}
	}

	Qhull qhull;
	qhull.runQhull("", dim, num_points, coords.data(), "Qt");

	// extract verts
	std::unordered_map<int, int> qhullIndexToLocalIndex;

	int localIndex = 0;
	for (auto v = qhull.vertexList().begin(); v != qhull.vertexList().end(); ++v) {
		QhullVertex vertex = *v;

		const double* pt = vertex.point().coordinates();

		Eigen::VectorXd p(dim);
		for (int d = 0; d < dim; ++d)
			p(d) = pt[d];
		
		poly.vertices.push_back(p);
		qhullIndexToLocalIndex[vertex.point().id()] = localIndex++;
	}

	//extract faces
	for (QhullFacet facet : qhull.facetList()) {
		if (!facet.isGood())
			continue;
		
		std::vector<int> face;

		QhullVertexSet vs = facet.vertices();
		for (auto vit = vs.begin(); vit != vs.end(); ++vit) {
			int qh_id = (*vit).point().id();
			face.push_back(qhullIndexToLocalIndex[qh_id]);
		}

		if (!face.empty())
			poly.faces.push_back(face);
	}
	return poly;
}

template <typename nodeType>
void betaPolytopes<nodeType>::computeHullForChart(Chart& chart, const std::vector<Eigen::VectorXd>& cloud_points) {
	std::vector<int> vertex_ids = collectChartVertexIndices(chart);
	chart.polytope = computeConvexHull(cloud_points, vertex_ids);
}

template <typename nodeType>
bool betaPolytopes<nodeType>::canMerge(Chart& A, Chart& B, const std::vector<Eigen::VectorXd>& cloud_points) {
	//collect vertices of the union of the charts
	std::unordered_set<int> union_vertex_set;

	auto addVertices = [&](Chart& C) {
		std::vector<int> vids = collectChartVertexIndices(C);
		for (int v : vids) union_vertex_set.insert(v);
	};

	addVertices(A);
	addVertices(B);

	std::vector<int> union_vertices(union_vertex_set.begin(), union_vertex_set.end());

	//compute hull of union
	Polytope union_hull = computeConvexHull(cloud_points, union_vertices);

	//convexity test (must contain exactly all union verts)
	if (union_hull.vertices.size() != union_vertex_set.size()) return false;

	return true;
}


template <typename nodeType>
void betaPolytopes<nodeType>::collect_ids(const std::vector<Strand>& strands, std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids, std::unordered_map<const Simplex*, int>& simplex_ids){
	int next_face_id = 0;
    int next_simplex_id = 0;

    for (const auto& strand : strands) {
        for (const auto& simplex : strand.simplices) {

            // Assign simplex ID if new
            if (!simplex_ids.count(simplex.get())) {
                simplex_ids[simplex.get()] = next_simplex_id++;
            }

            // Assign face IDs
            for (const auto& face : simplex->faces) {
				std::vector<unsigned> key = face -> vertices;
                if (!face_ids.count(key)) {
                    face_ids[key] = next_face_id++;
                }
            }
        }
    }
}

template <typename nodeType>
void betaPolytopes<nodeType>::write_faces_csv(const std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids)
	{
    std::ofstream file("../../python_tests/Polytopal_Development/faces.csv");
    file << "face_id,vertex_ids\n";

    for (const auto& [face, face_id] : face_ids) {
        file << face_id << ",\"";
        for (const auto& v : face) {
            file << v << " ";
        }
        file << "\"\n";
    }
}

template <typename nodeType>
void betaPolytopes<nodeType>::write_simplices_csv(const std::vector<Strand>& strands, std::unordered_map<std::vector<unsigned>, int, VectorHash>& face_ids, const std::unordered_map<const Simplex*, int>& simplex_ids){
	std::ofstream file("../../python_tests/Polytopal_Development/simplices.csv");
    file << "simplex_id,face_ids,strand_id\n";

    for (int strand_id = 0; strand_id < strands.size(); ++strand_id) {
        for (const auto& simplex : strands[strand_id].simplices) {

            int sid = simplex_ids.at(simplex.get());
            file << sid << ",\"";

            for (const auto& face : simplex->faces) {
				file << face_ids[face -> vertices] << " ";
            }

            file << "\"," << strand_id << "\n";
        }
    }
}

template <typename nodeType>
void betaPolytopes<nodeType>::export_strands_to_csv(const std::vector<Strand>& strands){
	std::unordered_map<std::vector<unsigned>, int, VectorHash> face_ids;
    std::unordered_map<const Simplex*, int> simplex_ids;

    collect_ids(strands, face_ids, simplex_ids);

    write_faces_csv(face_ids);
    write_simplices_csv(strands, face_ids, simplex_ids);
}

template <typename nodeType>
void betaPolytopes<nodeType>::exportAtlasStructure(const std::vector<Chart>& atlas) {
	std::cout << "Atlas contains " << atlas.size() << " charts\n\n";

	std::ofstream out("../../python_tests/Polytopal_Development/atlas.csv");

	if (!out.is_open()) {
		std::cerr << "Error: could not open file " << "atlas.csv" << "\n";
		return;
	}

	//csv header
	out << "chart_id,simplex_id,vertex_ids\n";

	size_t simplex_counter = 0;

	for (size_t c = 0; c < atlas.size(); ++c) {
		const Chart& chart = atlas[c];

		for (const auto& simplex : chart.simplices) {

			std::unordered_set<unsigned> vertex_ids;

			for (const auto& face : simplex -> faces) {
				for (unsigned vid : face -> vertices) {
					vertex_ids.insert(vid);
				}
			}

			out << c << "," << simplex_counter << ",";

			bool first = true;
			for (auto vid : vertex_ids) {
				if(!first) out << " ";
				out << vid;
				first = false;
			}

			out << "\n";

			simplex_counter++;
		}
	}

	out.close();
	std::cout << "atlas structure written\n";
}



template class betaPolytopes<simplexNode>;
template class betaPolytopes<alphaNode>;
template class betaPolytopes<witnessNode>;
