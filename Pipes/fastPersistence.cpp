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
#include "simplexTree.hpp"

unionFind::unionFind(int n) : rank(n, 0), parent(n, 0) {
	for(int i=0; i<n; i++) parent[i]=i;
}

int unionFind::find(int i){
	
	//
	if(i == parent[i]) return i;
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

//Hash the set by converting to a remapped index
long long fastPersistence::ripsIndex(std::set<unsigned>& simplex, binomialTable& bin){ 
	long long simplexIndex = 0;
	unsigned i = 0;
	auto it = simplex.begin();
	while(it != simplex.end()){
		simplexIndex += bin.binom(*it - shift, ++i);
		if(simplexIndex < 0) throw std::overflow_error("Binomial overflow");
		++it;
	}
	return simplexIndex;
}

unsigned fastPersistence::maxVertex(long long ripsIndex, unsigned high, unsigned low, unsigned k, binomialTable &bin){
	while(high > low){ //Binary search for the max vertex for this simplex
		unsigned mid = (high + low)/2;
		if(bin.binom(mid, k) <= ripsIndex) low = mid + 1;
		else high = mid;
	}
	return high - 1;
}

std::vector<unsigned> fastPersistence::getVertices(long long ripsIndex, int dim, unsigned n, binomialTable &bin){
	std::vector<unsigned> v;
	for(unsigned k = dim+1; k>0; k--){ //Get all vertices by repeated binary search for max vertex
		n = maxVertex(ripsIndex, n, k-1, k, bin);
		v.push_back(n);
		ripsIndex -= bin.binom(n, k);
	}
	return v;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//	FastPersistence: For computing the persistence pairs from simplicial complex:
//		1. See Bauer-19 for algorithm/description
pipePacket fastPersistence::runPipe(pipePacket inData){
	if(dim > 0)
		inData.complex->expandDimensions(dim + 1);	
	
	//Get all edges for the simplexArrayList or simplexTree
	std::vector<std::vector<simplexBase::simplexNode*>> edges = inData.complex->getAllEdges();

	if(edges.size() <= 1)
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
	unsigned nPts = inData.originalData.size();

	unionFind uf(nPts);
	shift = *edges[0][nPts-1]->simplex.begin();

	for(auto& edge : edges[1]){ //For each edge
		std::set<unsigned>::iterator it = edge->simplex.begin();
		
		//Find which connected component each vertex belongs to
		int v1 = uf.find(*it - shift), v2 = uf.find(*(++it) - shift); 
		
		//Edge connects two different components -> add to the MST
		if(v1 != v2){ 
			uf.join(v1, v2);
			pivots.push(edgeIndex);
			mstSize++;
			
			bettiBoundaryTableEntry des = { 0, 0, edge->weight, { *it - shift } };
			inData.bettiTable.push_back(des);
		}

		//Check if we've filled our MST and can break...
		if(mstSize >= edges[0].size()-1) break;

		edgeIndex++;
	}

	for(int i=0; i<inData.originalData.size(); i++){
		if(uf.find(i) == i){ //i is the name of a connected component
			//Each connected component has an open persistence interval
			bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {} };
			inData.bettiTable.push_back(des);
		}
	}
	
	//For higher dimensional persistence intervals
	//	
		//Build next dimension of ordered simplices, ignoring previous dimension pivots
		
		//Represent the ordered simplices as indexed sets; similar to the approach in indSimplexTree
		
		//Identify apparent pairs - i.e. d is the youngest face of d+1, and d+1 is the oldest coface of d
		//		This indicates a feature represents a persistence interval
		
		//Track V (reduction matrix) for each column j that has been reduced to identify the constituent 
		//		boundary simplices
		
	// if(true){
	// 	for(unsigned d = 1; d < dim; d++){
	// 		for(unsigned columnIndex = edges[d].size(); columnIndex-- != 0; ){
	// 			std::pair<std::set<unsigned>, double>& simplex = edges[d][columnIndex];
	// 			((simplexTree*)inData.complex)->getAllCofacets(simplex.first);
	// 		}
	// 	}
	// } else{

	binomialTable bin(nPts, dim+1);
	
	for(unsigned d = 1; d < dim && d < edges.size()-1; d++){
		
		//Convert our hashed set (ripsIndex) to the index from our complex
		std::unordered_map<long long, unsigned> indexConverter; 
		
		unsigned simplexIndex = 0;
		for(auto& simplex : edges[d+1]){
			indexConverter.insert(std::make_pair(ripsIndex(simplex->simplex, bin), simplexIndex));
			simplexIndex++;
		}

		std::unordered_map<unsigned, std::vector<long long>> v; 	//Store only the reduction matrix V and compute R implicity
		std::unordered_map<unsigned, unsigned> pivotPairs; 			//For each pivot, which column has that pivot
		std::priority_queue<unsigned> nextPivots; 					//Pivots for the next dimension

		//Iterate over columns to reduce in reverse order
		for(unsigned columnIndex = edges[d].size(); columnIndex-- != 0; ){ 
			
			simplexBase::simplexNode* simplex = edges[d][columnIndex];		//Pointer to the current simplex
			double simplexWeight = simplex->weight;										//Current simplex weight
			std::vector<unsigned> cofaceList; 											//Store cofaceList as a min heap
			std::vector<long long> columnV;												
			bool foundEmergentCandidate = false; 	//Found the lexicographically maximum cofacet with the same diameter

			//Not a pivot -> need to reduce
			if(pivots.top() != columnIndex){ 
				std::set<unsigned>::reverse_iterator it = simplex->simplex.rbegin();
				unsigned k = simplex->simplex.size() + 1;
				
				//Push back the remapped index for the current simplex
				long long index = ripsIndex(simplex->simplex, bin);
				columnV.push_back(index);

				//Try inserting other vertices into the simplex
				for(unsigned i=nPts; i-- != 0; ){ 
					
					if(it != simplex->simplex.rend() && i == *it - shift){ //Vertex i is already in the simplex
						//Now adding vertices less than i -> i is now the kth largest vertex in the simplex instead of the (k-1)th

						index -= bin.binom(i, k-1);
						index += bin.binom(i, k); //Recompute the index accordingly

						//Now need to check for the (k-1)th vertex in the simplex
						--k;
						
						//Check for the previous vertex in the simplex (it is a reverse iterator)
						++it;
					} else{
						
						
						auto cofacetIndex = indexConverter.find(index + bin.binom(i, k));
						if(cofacetIndex != indexConverter.end()){ //If this is a valid simplex, add it to the heap
							cofaceList.push_back(cofacetIndex->second);

							//If we haven't found an emergent candidate and the weight of the maximal cofacet is equal to the simplex's weight
							//		we have identified an emergent pair; at this point we can break because the interval is born and dies at the 
							//		same epsilon
							if(!foundEmergentCandidate && edges[d+1][cofacetIndex->second]->weight == simplex->weight){
								
								//Check to make sure the identified cofacet isn't a pivot
								if(pivotPairs.find(cofacetIndex->second) == pivotPairs.end()) 
									break; //Found an emergent cofacet pair -> we can break
									
								foundEmergentCandidate = true;
							}
						}
					}
				}

				//Never was able to iterate back to simplex.first.rend(); index holds the remapped coefficient
				//		build a heap using the coface list to reduce and store in V
				std::make_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());

				while(!cofaceList.empty()){
					
					unsigned pivotIndex; //Get minimum coface from heap
								
					
					while(!cofaceList.empty()){
						pivotIndex = cofaceList.front();
						
						//Rotate the heap
						std::pop_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
						
						
						cofaceList.pop_back();
						if(!cofaceList.empty() && pivotIndex == cofaceList.front()){ //Coface is in twice -> evaluates to 0 mod 2
							
							//Rotate the heap
							std::pop_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
							cofaceList.pop_back();
						} else{
							cofaceList.push_back(pivotIndex);
							
							std::push_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
							break;
						}
					}

					if(cofaceList.empty()){ //Column completely reduced
						break;
					} else if(pivotPairs.find(pivotIndex) == pivotPairs.end()){ //Column cannot be reduced
						nextPivots.push(pivotIndex);
						pivotPairs.insert({pivotIndex, columnIndex});

						// std::sort(columnV.begin(),columnV.end());
						// auto it = columnV.begin();
						// while(it != columnV.end()){
						// 	if((it+1)!=columnV.end() && *it==*(it+1)) ++it;
						// 	else{
						// 		v[columnIndex].push_back(*it);
						// 	}
						// 	++it;
						// }
						v[columnIndex] = columnV;

						if(edges[d][columnIndex]->weight != edges[d+1][pivotIndex]->weight){
							bettiBoundaryTableEntry des = { d, simplexWeight, edges[d+1][pivotIndex]->weight, std::set<unsigned>(cofaceList.begin(), cofaceList.end()) };
							inData.bettiTable.push_back(des);
						}

						break;
					} else{
						for(long long index : v[pivotPairs[pivotIndex]]){ //Reduce the column of R by computing the appropriate columns of D by enumerating cofacets
							columnV.push_back(index);
							std::vector<unsigned> facet = getVertices(index, d, nPts, bin); //Enumerate cofacets of this simplex like above
							std::vector<unsigned>::iterator it = facet.begin();
							unsigned k = facet.size() + 1;
							for(unsigned i = nPts; i-- != 0; ){ //Try inserting other vertices into the simplex
								if(it != facet.end() && i == *it){
									index -= bin.binom(i, k-1);
									index += bin.binom(i, k);
									--k;
									++it;
								} else{
									auto cofacetIndex = indexConverter.find(index + bin.binom(i, k));
									if(cofacetIndex != indexConverter.end()){ //If this is a valid simplex, add it to the heap
										cofaceList.push_back(cofacetIndex->second);
									}
								}
							}
						}
						std::make_heap(cofaceList.begin(), cofaceList.end(), std::greater<unsigned>());
					}
				}
				
				if(cofaceList.empty()){
					bettiBoundaryTableEntry des = { d, simplexWeight, maxEpsilon, {} };
					inData.bettiTable.push_back(des);
				}
			
			//Was a pivot, skip the evaluation and queue next pivot
			} else pivots.pop();
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
