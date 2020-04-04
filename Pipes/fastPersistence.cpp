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
#include <utility>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <queue>
#include "fastPersistence.hpp"

typedef std::priority_queue<unsigned, std::vector<unsigned>, std::greater<unsigned>> min_priority_queue;

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

binomialTable::binomialTable(unsigned n, unsigned k) : v(n+1, std::vector<long long>(k+1, 0)){ //Fast computation of binomials with precomputed table
	v[0][0] = 1;
	for(int i=1; i<=n; i++){
		v[i][0] = 1;
		for(int j=1; j<=k; j++){
			v[i][j] = v[i-1][j-1] + v[i-1][j]; //Pascal's Rule
			if(v[i][j] < 0) throw std::overflow_error("Binomial overflow");
		}
	}
}

long long binomialTable::binom(unsigned n, unsigned k){ //Return binomial coefficient
	if(k>n) return 0;
	return v[n][k];
}

// basePipe constructor
fastPersistence::fastPersistence(){
	pipeType = "FastPersistence";
	return;
}

long long fastPersistence::ripsIndex(std::set<unsigned>& simplex, binomialTable& bin){ //Hash the set by converting to the index described by the paper
	long long simplexIndex = 0;
	unsigned i = 0;
	auto it = simplex.begin();
	while(it != simplex.end()){
		simplexIndex += bin.binom(*it, ++i);
		if(simplexIndex < 0) throw std::overflow_error("Binomial overflow");
		++it;
	}
	return simplexIndex;
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
	binomialTable bin(inData.originalData.size(), dim+1);

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
		std::unordered_map<long long, unsigned> indexConverter; //Convert our hashed set (ripsIndex) to the index from our complex
		unsigned simplexIndex = 0;
		for(auto simplex : edges[d+1]){
			indexConverter.insert(std::make_pair(ripsIndex(simplex.first, bin), simplexIndex));
			simplexIndex++;
		}

		std::unordered_map<unsigned, std::vector<unsigned>> boundary;
		std::unordered_map<unsigned, unsigned> pivotPairs; //For each pivot, which column has that pivot
		std::priority_queue<unsigned> nextPivots; //Pivots for the next dimension

		for(unsigned columnIndex = edges[d].size()-1; columnIndex-- != 0; ){ //Iterate over columns to reduce in reverse order
			std::pair<std::set<unsigned>, double> simplex = edges[d][columnIndex];
			double simplexWeight = simplex.second;
			std::vector<unsigned> cofaceList; //Store cofaceList as a min heap

			if(pivots.top() != columnIndex){ //Not a pivot -> need to reduce
				unsigned cofacetIndex = 0;

				std::set<unsigned> cofacet;
				for(unsigned i=0; i<inData.originalData.size(); i++){ //Try inserting other vertices into the simplex
					cofacet = simplex.first;
					if(cofacet.insert(i).second){ //New vertex was added to the simplex
						auto simplexIndex = indexConverter.find(ripsIndex(cofacet, bin)); //Convert our set to its index using the ripsIndex as an intermediate hash
						if(simplexIndex != indexConverter.end()){ //If this is a valid simplex, add it to the heap
							cofaceList.push_back(simplexIndex->second);
							std::push_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
						}
					}
				}

				if(cofaceList.empty()){
					bettis += std::to_string(d) + "," + std::to_string(simplexWeight) +"," + std::to_string(maxEpsilon) + "\n";
					bettiBoundaryTableEntry des = { d, simplexWeight, maxEpsilon, {} };
					inData.bettiTable.push_back(des);
				}

				while(!cofaceList.empty()){
					unsigned pivotIndex; //Get minimum coface from heap
					while(true){
						pivotIndex = cofaceList.front();
						std::pop_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
						cofaceList.pop_back();
						if(pivotIndex == cofaceList.front()){ //Coface is in twice -> evaluates to 0 mod 2
							std::pop_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
							cofaceList.pop_back();
						} else break;
					}

					if(pivotPairs.find(pivotIndex) == pivotPairs.end()){
						nextPivots.push(pivotIndex);
						pivotPairs.insert({pivotIndex, columnIndex});	
						boundary[columnIndex] = cofaceList;
						
						if(edges[d][columnIndex].second != edges[d+1][pivotIndex].second){
							bettis += std::to_string(d) + "," + std::to_string(simplexWeight) +"," + std::to_string(edges[d+1][pivotIndex].second) + "\n";
							bettiBoundaryTableEntry des = { d, simplexWeight, edges[d+1][pivotIndex].second, std::set<unsigned>(cofaceList.begin(), cofaceList.end()) };
							inData.bettiTable.push_back(des);
						}

						break;
					} else {
						cofaceList.insert(cofaceList.end(), boundary[pivotPairs[pivotIndex]].begin(), boundary[pivotPairs[pivotIndex]].end()); //Lazy set XOR
						std::make_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
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

