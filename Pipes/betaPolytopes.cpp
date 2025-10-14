/*
 * betaSubSkeletonComplexPipe hpp + cpp extend the basePipe class for calculating the
 * beta Skeleton Based Complex generation for data input
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
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
	std::vector<Simplex> mesh_structs; //copy in here to have struct features in each simplex
	//mesh_structs.reserve(this -> dim + 1); potential optimization
	std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq> facelist; //stores <face, #of incident simplex>
	std::vector<Strand> strands;

	for(auto x:dsimplexmesh){
	  for(auto y:x)
		std::cout<<y<<" ";
	std::cout<<std::endl;
	}
	
	std::cout<<"We will generate Polytopes here"<<std::endl;
	
	for(const auto& x:dsimplexmesh) {
		auto current_simplex = std::make_shared<Simplex>(this->dim);
		for(int i = 0; i < x.size(); ++i) {
			//std::vector<unsigned> current_face;
			//current_face.reserve(x.size() - 1);
			auto current_face = std::make_shared<Face>();
			for(int j = 0; j < x.size(); ++j) {
				if (i != j){
					current_face -> verticies.push_back(x[j]);
				}
			}
			current_simplex -> faces.push_back(current_face);
			current_face -> adjacent_simplicies.push_back(current_simplex);
			auto [it, inserted] = facelist.emplace(current_face, 1);
			if(!inserted) {
				it -> second++;
			}
		}
		mesh_structs.push_back(*current_simplex);
	}

	for(auto simplex:mesh_structs) {
		if(simplex.visited == false) {
			Strand strand;
			flood_fill(strand, simplex, facelist);
			strands.push_back(strand);
		}
	}

//TESTS

	std::cout << "print strands\n";
	for(const auto& strand:strands) {
		for(const auto& simplex:strand.simplicies) {
			for(const auto& face:simplex.faces) {
				for(const auto& vert:face -> verticies) {
					std::cout << vert << " ";
				}
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

/*
	std::cout << "print simplicies as face lists\n";
	for(const auto& simplex : mesh_structs) {
		for(auto face:simplex.faces) {
			for(auto x:face){
				std::cout << x << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}



	std::cout << "dump unordered_map" << std::endl;
	std::cout << "facelist size: " << facelist.size() << std::endl;
	for (auto& [f, deg] : facelist) {
		std::cout << "Face { ";
		for (auto v : f->verticies)
			std::cout << v << " ";
		std::cout << "} degree=" << deg << std::endl;
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
void betaPolytopes<nodeType>::flood_fill(Strand& strand, Simplex& simplex, const std::unordered_map<std::shared_ptr<Face>, unsigned, FacePtrHash, FacePtrEq>& facelist) {
	if(simplex.visited) return;
	strand.simplicies.push_back(simplex);
	simplex.visited = true;
	for(const auto& face:simplex.faces) {
		if(facelist.at(face)==this -> dim) {
			for(auto& current_simplex:face -> adjacent_simplicies) {
				flood_fill(strand, *current_simplex, facelist);
			}
		}
	}
	return;
}



template class betaPolytopes<simplexNode>;
template class betaPolytopes<alphaNode>;
template class betaPolytopes<witnessNode>;
