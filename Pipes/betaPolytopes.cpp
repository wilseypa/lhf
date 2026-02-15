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
/*
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

	for(auto simplex:mesh_structs) {
		if(simplex -> visited == false) {
			Strand strand; //whenever last recursion ends (last strand fully enumerated), find a new unvisited simplex to enumerate new strand
			flood_fill(strand, simplex, facelist);
			strands.push_back(strand);
		}
	}

	export_strands_to_csv(strands);

	std::cout << "Generating atlases from strands..." << std::endl;

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
		generateAtlasForStrand(strand, cloud_points, 2, 0.001, facelist);
	}

	//assigning ids and exporting atlases to files
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

		for (size_t chart_id = 0; chart_id < strand.atlas.size(); ++chart_id) {
			const auto& chart = strand.atlas[chart_id];

			for (const auto& simplex : chart.simplices) {
				atlasFile
					<< strand_id << ","
					<< chart_id << ","
					<< simplex_ids[simplex.get()]
					<< "\n";
			}
		}
	}

	atlasFile.close();

	std::cout << "Atlas CSV written successfully." << std::endl;


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
void betaPolytopes<nodeType>::generateAtlasForStrand(Strand& strand, const std::vector<Eigen::VectorXd>& cloud_points, int intrinsic_dim, double distortion_threshold, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> &facelist) {
	//reset simplices visited
	for(auto& simplex: strand.simplices) {
		simplex->visited = false;
	}

	//iterate through until all simplices in strand are assigned
	for (auto& seed : strand.simplices) {
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

			//get neighbors
			for(const auto& face:current -> faces) {
				if(facelist.at(face) != 2) continue;

				std::unordered_set<std::shared_ptr<Simplex>> unique_neighbors(
					face->adjacent_simplices.begin(),
					face->adjacent_simplices.end()
				);

				for(auto& neighbor : unique_neighbors) {
					if(neighbor->visited) continue;

					if(chart.tryAddSimplex(neighbor, cloud_points)) {
						neighbor->visited = true;
						q.push(neighbor);
					}
				}
			}
		}
		//store chart in atlas
		strand.atlas.push_back(chart);
	}
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



template class betaPolytopes<simplexNode>;
template class betaPolytopes<alphaNode>;
template class betaPolytopes<witnessNode>;
