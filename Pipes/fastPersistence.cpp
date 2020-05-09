/*
 * persistencePairs hpp + cpp extend the basePipe class for calculating the
 * persistence pairs numbers from a complex
 *
 */

#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <queue>
#include "fastPersistence.hpp"
#include "utils.hpp"

unionFind::unionFind(int n) : rank(n, 0), parent(n, 0) {
	for(int i=0; i<n; i++) parent[i]=i;
}

int unionFind::find(int i){
	if(i == parent[i]) return i; //Found name of the component
	parent[i] = find(parent[i]); //Path Compression
	return parent[i];
}

bool unionFind::join(int x, int y){ //Union by rank
	x = find(x);
	y = find(y);
	if(x == y) return false;
	if(rank[x] == rank[y]){
		rank[y]++;
		parent[x] = y;
	} else if(rank[x] < rank[y]){
		parent[x] = y;
	} else{
		parent[y] = x;
	}
	return true;
}

// basePipe constructor
fastPersistence::fastPersistence(){
	pipeType = "FastPersistence";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//	FastPersistence: For computing the persistence pairs from simplicial complex:
//		1. See Bauer-19 for algorithm/description
pipePacket fastPersistence::runPipe(pipePacket inData){
	if(dim > 0) inData.complex->expandDimensions(dim + 1);

	//Get all edges for the simplexArrayList or simplexTree
	std::vector<std::set<simplexNode*, cmpByWeight>> edges = inData.complex->getAllEdges();


	if(edges.size() <= 1) return inData;

	//Some notes on fast persistence:

	//	-Vectors need to be stored in a lexicograhically ordered set of decreasing (d+1)-tuples (e.g. {3, 1, 0})
	//		-These vectors are replaced with their indices (e.g. {{2,1,0} = 0, {3,1,0} = 1, {3,2,0} = 2, etc.})

	//	-Boundary matrix stores reduced collection of cofaces for each column (indexed)

	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();

	//Get all dim 0 persistence intervals
	//Kruskal's minimum spanning tree algorithm
	//		Track the current connected components (uf)
	//		Check edges; if contained in a component, ignore
	//			if joins components, add them
	//		Until all edges evaluated or MST found (size - 1)

	std::vector<simplexNode*> pivots; //Store identified pivots
	unsigned mstSize = 0;
	unsigned nPts = inData.originalData.size();

	unionFind uf(nPts);
	//shift = *(*edges[0].begin())->simplex.begin();
	shift = 0;

	for(auto edgeIter = edges[1].begin(); edgeIter != edges[1].end(); edgeIter++){
		std::set<unsigned>::iterator it = (*edgeIter)->simplex.begin();

		//Find which connected component each vertex belongs to
		int v1 = uf.find(*it - shift),  v2 = uf.find(*(++it) - shift);

		//Edge connects two different components -> add to the MST
		if(v1 != v2){
			uf.join(v1, v2);
			mstSize++;

			simplexNode* temp = new simplexNode;
			temp->simplex = (*edgeIter)->simplex;
			temp->weight = (*edgeIter)->weight;
			pivots.push_back(temp);

			bettiBoundaryTableEntry des = { 0, 0, (*edgeIter)->weight, {}, {temp} };
			inData.bettiTable.push_back(des);
		}

		//Check if we've filled our MST and can break
		if(mstSize >= edges[0].size()-1) break;
	}

	for(int i=0; i<inData.originalData.size(); i++){
		if(uf.find(i) == i){ //i is the name of a connected component
			//Each connected component has an open persistence interval
			bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {}, {} };
			inData.bettiTable.push_back(des);
		}
	}


	//For higher dimensional persistence intervals
	//
		//Build next dimension of ordered simplices, ignoring previous dimension pivots

		//Represent the ordered simplices as indexed sets

		//Identify apparent pairs - i.e. d is the youngest face of d+1, and d+1 is the oldest coface of d
		//		This indicates a feature represents a trivial persistence interval

		//Track V (reduction matrix) for each column j that has been reduced to identify the constituent
		//		boundary simplices

	for(unsigned d = 1; d < dim && d < edges.size()-1; d++){
		std::sort(pivots.begin(), pivots.end(), cmpByWeight());
		std::vector<simplexNode*>::iterator it = pivots.begin();

		//Convert our hashed set (ripsIndex) to the index from our complex
		//std::unordered_map<long long, unsigned> indexConverter;

		//unsigned simplexIndex = 0;
		//for(auto& simplex : edges[d+1]){
			//indexConverter.insert(std::make_pair(ripsIndex(simplex->simplex, bin), simplexIndex));
			//simplexIndex++;
		//}

		std::vector<simplexNode*> nextPivots;	 					//Pivots for the next dimension
		std::unordered_map<unsigned, std::vector<unsigned>> v;				//Store only the reduction matrix V and compute R implicity
		std::unordered_map<simplexNode*, unsigned> pivotPairs;	//For each pivot, which column has that pivot

		unsigned columnIndex = edges[d].size() - 1; //Keep track of the column index for the pivots (for now)

		//Iterate over columns to reduce in reverse order
		for(auto columnIndexIter = edges[d].rbegin(); columnIndexIter != edges[d].rend(); columnIndexIter++, columnIndex--){
			simplexNode* simplex = (*columnIndexIter);		//The current simplex

			//Not a pivot -> need to reduce
			if((*it)->weight != simplex->weight && (*it)->simplex != simplex->simplex){
				//Get all cofacets using emergent pair optimization
				std::vector<simplexNode*> cofaceList = inData.complex->getAllCofacets(simplex->simplex, simplex->weight, pivotPairs);
				std::vector<unsigned> columnV;	//Reduction column of matrix V
				columnV.push_back(columnIndex); //Initially V=I -> 1's along diagonal

				//Build a heap using the coface list to reduce and store in V
				std::make_heap(cofaceList.begin(), cofaceList.end(), cmpBySecond());

				while(true){
					simplexNode* pivot;
					while(!cofaceList.empty()){
						pivot = cofaceList.front();

						//Rotate the heap
						std::pop_heap(cofaceList.begin(), cofaceList.end(), cmpBySecond());
						cofaceList.pop_back();

						if(!cofaceList.empty() && pivot == cofaceList.front()){ //Coface is in twice -> evaluates to 0 mod 2

							//Rotate the heap
							std::pop_heap(cofaceList.begin(), cofaceList.end(), cmpBySecond());
							cofaceList.pop_back();
						} else{

							cofaceList.push_back(pivot);
							std::push_heap(cofaceList.begin(), cofaceList.end(), cmpBySecond());
							break;
						}
					}

					if(cofaceList.empty()){ //Column completely reduced
						break;
					} else if(pivotPairs.find(pivot) == pivotPairs.end()){ //Column cannot be reduced
						pivotPairs.insert({pivot, columnIndex});
						nextPivots.push_back(pivot);

						// std::sort(columnV.begin(), columnV.end());
						// auto it = columnV.begin();
						// while(it != columnV.end()){
						// 	if((it+1) != columnV.end() && *it==*(it+1)) ++it;
						// 	else v[columnIndex].push_back(*it);
						// 	++it;
						// }

						v[columnIndex] = columnV;

						if(simplex->weight != pivot->weight){
							bettiBoundaryTableEntry des = { d, simplex->weight, pivot->weight, {}, cofaceList };

							inData.bettiTable.push_back(des);
						}

						break;
					} else{ //Reduce the column of R by computing the appropriate columns of D by enumerating cofacets
						for(unsigned i : v[pivotPairs[pivot]]){
							auto temp = edges[d].begin();
							std::advance(temp, i);

							std::vector<simplexNode*> cofaces = inData.complex->getAllCofacets((*temp)->simplex);
							cofaceList.insert(cofaceList.end(), cofaces.begin(), cofaces.end());
						}
						std::make_heap(cofaceList.begin(), cofaceList.end(), cmpBySecond());
					}
				}

			//Was a pivot, skip the evaluation and queue next pivot
			} else ++it;
		}

		pivots = nextPivots;
	}

	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();

	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;

	//Output the time and memory used for this pipeline segment
	ut.writeDebug("persistence","Bettis executed in " + std::to_string(elapsed.count()/1000.0) + " seconds (physical time)");;

	return inData;
}



// outputData -> used for tracking each stage of the pipeline's data output without runtime
void fastPersistence::outputData(pipePacket inData){
	std::ofstream file;
	if(fnmod.size() > 0)
		file.open("output/"+pipeType+"_bettis_output"+fnmod+".csv");
	else
		file.open("output/" + pipeType + "_bettis_output.csv");

	file << inData.bettiOutput;

	file.close();

	file.open("output/tArray.csv");

	file << "Dim,Birth,Death,Simplex\n";
	for(auto tStruct : inData.bettiTable){
		file << tStruct.bettiDim << "," << tStruct.birth << "," << tStruct.death << ",";
		for(auto index : tStruct.boundaryPoints)
			file << index << " ";
		file << "\n";
	}

	file.close();

	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool fastPersistence::configPipe(std::map<std::string, std::string> configMap){
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = configMap["outputFile"].c_str();

	ut = utils(strDebug, outputFile);

	pipe = configMap.find("dimensions");
	if(pipe != configMap.end())
		dim = std::atoi(configMap["dimensions"].c_str());
	else return false;

	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return false;

	pipe = configMap.find("fn");
	if(pipe != configMap.end())
		fnmod = configMap["fn"];

	configured = true;
	ut.writeDebug("fastPersistence","Configured with parameters { dim: " + configMap["dimensions"] + ", complexType: " + configMap["complexType"] + ", eps: " + configMap["epsilon"]);
	ut.writeDebug("fastPersistence","\t\t\t\tdebug: " + strDebug + ", outputFile: " + outputFile + " }");

	return true;
}
