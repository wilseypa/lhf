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

// basePipe constructor
template<typename T>
fastPersistence<T>::fastPersistence(){
	this->pipeType = "FastPersistence";
	return;
}

template <class T>
template <class simplexNodePointer, class comp>
std::vector<simplexNodePointer> fastPersistence<T>::persistenceByDimension(pipePacket<T>& inData, std::vector<simplexNodePointer> edges, std::vector<simplexNodePointer> pivots, unsigned dimension, comp compStruct, std::string mode, bool recordIntervals){
	std::sort(edges.begin(), edges.end(), compStruct);
	std::sort(pivots.begin(), pivots.end(), compStruct);

	typename std::vector<simplexNodePointer>::iterator it = pivots.begin();

	std::vector<simplexNodePointer> nextPivots;	 	//Pivots for the next dimension
	std::unordered_map<simplexNodePointer, std::vector<simplexNodePointer>> v;	//Store only the reduction matrix V and compute R implicity
	std::unordered_map<simplexNodePointer, simplexNodePointer> pivotPairs;	//For each pivot, which column has that pivot
	
	//Iterate over columns to reduce in reverse order
	for(auto columnIndexIter = edges.begin(); columnIndexIter != edges.end(); columnIndexIter++){
		simplexNodePointer simplex = (*columnIndexIter);	//The current simplex

		//Not a pivot -> need to reduce
		if(it == pivots.end() || (*it)->weight != simplex->weight || (*it)->simplex != simplex->simplex){
		//	std::cout<<mode<<" "<<simplicialComplex<<" "<<complexType<<std::endl;
			//Get all cofacets using emergent pair optimization
			std::vector<simplexNodePointer> faceList = (mode == "homology" ? inData.complex->getAllFacets_P(simplex) : ((this->simplicialComplex == "alpha" && this->complexType == "simplexArrayList")? inData.complex->getAllDelaunayCofacets(simplex):inData.complex->getAllCofacets(simplex->simplex, simplex->weight, pivotPairs, true)));
		
			std::vector<simplexNodePointer> columnV;	//Reduction column of matrix V
			columnV.push_back(simplex); //Initially V=I -> 1's along diagonal

			//Build a heap using the coface list to reduce and store in V
			std::make_heap(faceList.begin(), faceList.end(), compStruct);

			while(true){
				simplexNodePointer pivot;

				while(!faceList.empty()){
					pivot = faceList.front();

					//Rotate the heap
					std::pop_heap(faceList.begin(), faceList.end(), compStruct);
					faceList.pop_back();

					if(!faceList.empty() && pivot->hash == faceList.front()->hash){ //Coface is in twice -> evaluates to 0 mod 2

						//Rotate the heap
						std::pop_heap(faceList.begin(), faceList.end(), compStruct);
						faceList.pop_back();
					} else{

						faceList.push_back(pivot);
						std::push_heap(faceList.begin(), faceList.end(), compStruct);
						break;
					}
				}

				if(faceList.empty()){ //Column completely reduced
					break;
				} else if(pivotPairs.find(pivot) == pivotPairs.end()){ //Column cannot be reduced
					pivotPairs.insert({pivot, simplex});
					nextPivots.push_back(pivot);
					
					std::sort(columnV.begin(), columnV.end());
					auto it = columnV.begin();
					while(it != columnV.end()){
						if((it+1) != columnV.end() && (*it)==*(it+1)) ++it;
						else v[simplex].push_back(*it);
						++it;
					}

					if(recordIntervals && simplex->weight != pivot->weight){
						bettiBoundaryTableEntry des = { dimension, std::min(pivot->weight, simplex->weight), std::max(pivot->weight, simplex->weight), this->ut.extractBoundaryPoints(v[simplex]) };
						inData.bettiTable.push_back(des);
					}

					break;
				} else{ 
					//Reduce the column of R by computing the appropriate columns of D by enumerating cofacets
					for(simplexNodePointer simp : v[pivotPairs[pivot]]){
						columnV.push_back(simp);
						std::vector<simplexNodePointer> faces = (mode == "homology" ? inData.complex->getAllFacets_P(simp) : ((this->simplicialComplex == "alpha" && this->complexType =="simplexArrayList")? inData.complex->getAllDelaunayCofacets(simp):inData.complex->getAllCofacets(simp->simplex)));
					
						faceList.insert(faceList.end(), faces.begin(), faces.end());
					}
					std::make_heap(faceList.begin(), faceList.end(), compStruct);
				}
			}
		//Was a pivot, skip the evaluation and queue next pivot
		} else ++it;
	}

	return nextPivots;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//	FastPersistence: For computing the persistence pairs from simplicial complex:
//		1. See Bauer-19 for algorithm/description
template <class T>
void fastPersistence<T>::runPipe(pipePacket<T> &inData){
	//Get all edges for the simplexArrayList or simplexTree
	std::vector<std::set<std::shared_ptr<T>, cmpByWeight<std::shared_ptr<T>>>> edges = inData.complex->getAllEdges();

	if(edges.size() <= 1) return;

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

	//For streaming data, indices will not be 0-N; instead sparse
	//	So in streaming, create a hash map to quickly lookup points

	std::unordered_map<unsigned, unsigned> mappedIndices;	//Store a map of the indices for MST
	std::vector<std::shared_ptr<T>> pivots; //Store identified pivots
	unsigned mstSize = 0;
	unsigned nPts = inData.workData.size();

	unionFind uf(nPts);

	for(auto edgeIter = edges[1].begin(); edgeIter != edges[1].end(); edgeIter++){
		std::set<unsigned>::iterator it = (*edgeIter)->simplex.begin();

		//Find which connected component each vertex belongs to
		//	Use a hash map to track insertions for streaming or sparse indices
		if( mappedIndices.size() == 0 || mappedIndices.find(*it) == mappedIndices.end() ) mappedIndices.insert( std::make_pair(*it, mappedIndices.size()) );
		int c1 = uf.find(mappedIndices.find(*it)->second);
		it++;
		if( mappedIndices.find(*it) == mappedIndices.end() ) mappedIndices.insert( std::make_pair(*it, mappedIndices.size()) );
		int c2 = uf.find(mappedIndices.find(*it)->second);

		//Edge connects two different components -> add to the MST
		if(c1 != c2){
			uf.join(c1, c2);
			mstSize++;

			pivots.push_back((*edgeIter));

			bettiBoundaryTableEntry des = { 0, 0, (*edgeIter)->weight, (*edgeIter)->simplex };
			inData.bettiTable.push_back(des);
		}

		//Check if we've filled our MST and can break
		if(mstSize >= edges[0].size()-1) break;
	}

	// std::cout << "mappedIndices.size = " << mappedIndices.size() << '\n';

	for(int i=0; i<inData.workData.size(); i++){
		if(uf.find(i) == i){ //i is the name of a connected component
			//Each connected component has an open persistence interval
			bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {} };
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

	bool involuted = (inv == "true");

	for(unsigned d = 1; d < dim && d < edges.size()-1; d++){
		inData.complex->prepareCofacets(d);

		pivots = persistenceByDimension(inData, std::vector<std::shared_ptr<T>>(edges[d].begin(), edges[d].end()), pivots, d, sortReverseLexicographic(), "cohomology", !involuted);

		//To recover the representative cycles from the cocycles, we compute homology on just the pivot columns
		if(involuted){
			inData.complex->prepareFacets(d);
			persistenceByDimension(inData, pivots, std::vector<std::shared_ptr<T>>(), d, sortLexicographic(), "homology", true);
		}
	}

	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();

	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;

	//Output the time and memory used for this pipeline segment
	this->ut.writeDebug("persistence","Bettis executed in " + std::to_string(elapsed.count()/1000.0) + " seconds (physical time)");;

	return;
}



// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <class T>
void fastPersistence<T>::outputData(pipePacket<T> &inData){
	std::ofstream file;
	if(this->fnmod.size() > 0)
		file.open("output/"+this->pipeType+"_bettis_output"+this->fnmod+".csv");
	else
		file.open("output/" + this->pipeType + "_bettis_output.csv");

	for(auto row : inData.bettiTable)
		file << std::to_string(row.bettiDim) << "," << std::to_string(row.birth) << "," << std::to_string(row.death) << std::endl;
	
	file.close();

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
template <class T>
bool fastPersistence<T>::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();

	this->ut = utils(strDebug, this->outputFile);

	pipe = configMap.find("involuted");
	if(pipe != configMap.end())
		this->inv = configMap["involuted"];

	pipe = configMap.find("dimensions");
	if(pipe != configMap.end())
		this->dim = std::atoi(configMap["dimensions"].c_str());
	else return false;

	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		this->maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return false;

	pipe = configMap.find("fn");
	if(pipe != configMap.end())
		this->fnmod = configMap["fn"];

	pipe = configMap.find("simplicialComplex");
	if(pipe != configMap.end())
		this->simplicialComplex = configMap["simplicialComplex"];

	pipe = configMap.find("complexType");
	if(pipe != configMap.end())
		this->complexType = configMap["complexType"];
	this->configured = true;
	this->ut.writeDebug("fastPersistence","Configured with parameters { dim: " + configMap["dimensions"] + ", complexType: " + configMap["complexType"] + ", eps: " + configMap["epsilon"]);
	this->ut.writeDebug("fastPersistence","\t\t\t\tdebug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}
