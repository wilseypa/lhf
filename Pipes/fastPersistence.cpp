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
#include <cmath>
#include <list>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <queue>
#include "fastPersistence.hpp"

unionFind::unionFind(int n) : rank(n, 0), parent(n, 0) {
	for(int i=0; i<n; i++) parent[i]=i;
}

int unionFind::find(int i){
	if(i == parent[i]) return i;
	parent[i] = find(parent[i]); //Path Compression
	return parent[i];
	// while(parent[i] != i){
	// 	parent[i] = parent[parent[i]];
	// 	i = parent[i];
	// }
	// return i;
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
	if(dim > 0)
		inData.complex->expandDimensions(dim + 1);	
	
	std::string bettis = "";
	
	//Get all edges for the simplexArrayList or simplexTree
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges = inData.complex->getAllEdges(maxEpsilon);

	if(edges.size() == 0)
		return inData;
		
	//Some notes on fast persistence:
	
	//	-Vectors need to be stored in a lexicograhically ordered set of decreasing (d+1)-tuples (e.g. {3, 1, 0})
	//		-These vectors are replaced with their indices (e.g. {{2,1,0} = 0, {3,1,0} = 1, {3,2,0} = 2, etc.})
	
	//	-Boundary matrix stores reduced collection of cofaces for each column (indexed)
	
	
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	//Get all dim 0 persistence intervals
	//Kruskal's minimum spanning tree algorithm
	//		Track the current connected components (uf)
	//		Check edges; if contained in current set ignore
	//			if joins multiple sets, join all sets
	//		Until all edges evaluated or MST found (size - 1)	
	
	std::priority_queue<unsigned> pivots; //Store identified pivots
	unsigned edgeIndex = 0;
	unsigned mstSize = 0;

	unionFind uf(inData.originalData.size());

	for(auto edge : edges[1]){ //For each edge
		std::set<unsigned>::iterator it = edge.first.begin();
		int v1 = uf.find(*it), v2 = uf.find(*(++it)); //Find which connected component each vertex belongs to
		if(v1 != v2){ //Edge connects two different components -> add to the MST
			uf.join(v1, v2);
			pivots.push(edgeIndex);
			mstSize++;
			bettis += "0,0," + std::to_string(edge.second) + "\n";
			bettiBoundaryTableEntry des = { 0, 0, edge.second, { *it } };
			inData.bettiTable.push_back(des);
		}

		//Check if we've filled our MST and can break...
		if(mstSize == edges[0].size()-1) break;

		edgeIndex++;
	}
	
	bettis += "0,0," + std::to_string(maxEpsilon) + "\n";
	bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {} };
	inData.bettiTable.push_back(des);
	
	//For higher dimensional persistence intervals
	//	
		//Build next dimension of ordered simplices, ignoring previous dimension pivots
		
		//Represent the ordered simplices as indexed sets; similar to the approach in indSimplexTree
		
		//Identify apparent pairs - i.e. d is the youngest face of d+1, and d+1 is the oldest coface of d
		//		This indicates a feature represents a persistence interval
		
		//Track V (reduction matrix) for each column j that has been reduced to identify the constituent 
		//		boundary simplices; we may not track this currently but will eventually
		
		
	for(unsigned d = 1; d < dim && d < edges.size()-1; d++){
		std::unordered_map<unsigned, std::set<unsigned>> boundary;
		std::unordered_map<unsigned, unsigned> pivotPairs; //For each pivot, which column has that pivot
		std::priority_queue<unsigned> nextPivots; //Pivots for the next dimension

		for(unsigned columnIndex = edges[d].size()-1; columnIndex-- != 0; ){ //Iterate over columns to reduce in reverse order
			std::pair<std::set<unsigned>, double> simplex = edges[d][columnIndex];
			std::set<unsigned> cofaceList;

			if(pivots.top() != columnIndex){ //Not a pivot -> need to reduce
				unsigned cofacetIndex = 0;

				// for(unsigned cofacetIndex = edges[d+1].size()-1; cofacetIndex-- != 0; ){
					// auto cofacetCandidate = edges[d+1][cofacetIndex];
				for(std::pair<std::set<unsigned>,double> cofacetCandidate : edges[d+1]){ //Iterate over possible cofacets
					bool isCofacet = std::includes(cofacetCandidate.first.begin(), cofacetCandidate.first.end(), simplex.first.begin(), simplex.first.end());

					if(isCofacet){
						cofaceList.insert(cofacetIndex);
						// if(check_apparent && cofacetCandidate.second == simplex.second){
						// 	if(pivotPairs.find(cofacetIndex) == pivotPairs.end()) break;
						// 	check_apparent = false;
						// // 	break; //Apparent Pair -> Don't need to enumerate the rest of the cofacets
						// }
					}
					++cofacetIndex;
				}

				if(cofaceList.empty()){
					bettis += std::to_string(d) + "," + std::to_string(simplex.second) +"," + std::to_string(maxEpsilon) + "\n";
					bettiBoundaryTableEntry des = { d, simplex.second, maxEpsilon, {} };
					inData.bettiTable.push_back(des);
				}

				while(!cofaceList.empty()){
					unsigned pivotIndex = *cofaceList.begin();

					if(pivotPairs.find(pivotIndex) == pivotPairs.end()){
						nextPivots.push(pivotIndex);
						pivotPairs.insert({pivotIndex, columnIndex});	
						boundary[columnIndex] = cofaceList;
						
						if(edges[d][columnIndex].second != edges[d+1][pivotIndex].second){
							bettis += std::to_string(d) + "," + std::to_string(simplex.second) +"," + std::to_string(edges[d+1][pivotIndex].second) + "\n";
							bettiBoundaryTableEntry des = { d, simplex.second, edges[d+1][pivotIndex].second, cofaceList };
							inData.bettiTable.push_back(des);
						}

						break;
					} else {
						cofaceList = ut.setXOR(cofaceList, boundary[pivotPairs[pivotIndex]]);
					}
				}

			} else pivots.pop();
		}
		
		pivots = nextPivots;
	}
			
	std::cout << bettis << std::endl;
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();
	
	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
	
	//Output the time and memory used for this pipeline segment
	ut.writeDebug("persistence","Bettis executed in " + std::to_string(elapsed.count()/1000.0) + " seconds (physical time)");;
	
	//Print the bettis
	//if(debug){
		//std::cout << std::endl << bettis << std::endl;
		//for(auto a : inData.bettiTable){
			//std::cout << a.bettiDim << "," << a.birth << "," << a.death << ",";
			//ut.print1DVector(a.boundaryPoints);
		//}
	//}
		
	inData.bettiOutput = bettis;
		
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

