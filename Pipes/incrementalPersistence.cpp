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
#include <exception>
#include <unordered_map>
#include <queue>
#include "incrementalPersistence.hpp"
#include "simplexArrayList.hpp"
#include "utils.hpp"

// basePipe constructor
template<typename nodeType>
incrementalPersistence<nodeType>::incrementalPersistence(){
	this->pipeType = "IncrementalPersistence";
	return;
}

template <class nodeType>
template <typename simplexNodePointer, typename comp>
std::vector<simplexNodePointer> incrementalPersistence<nodeType>::incrementalByDimension(pipePacket<nodeType>& inData, std::vector<simplexNodePointer>& edges, std::vector<simplexNodePointer> pivots, unsigned dimension, comp compStruct, std::string mode, bool recordIntervals){
	std::sort(edges.begin(), edges.end(), compStruct);
	std::sort(pivots.begin(), pivots.end(), compStruct);

	typename std::vector<simplexNodePointer>::iterator it = pivots.begin();

	std::vector<simplexNodePointer> nextPivots;	 	//Pivots for the next dimension
	std::unordered_map<simplexNodePointer, std::vector<simplexNodePointer>> v;	//Store only the reduction matrix V and compute R implicity
	std::unordered_map<long long, simplexNodePointer> pivotPairs;	//For each pivot, which column has that pivot

	simplexArrayList<nodeType>* complex;
	if(inData.complex->simplexType == "simplexArrayList" || inData.complex->simplexType == "alphaComplex"){
		complex = (simplexArrayList<nodeType>*) inData.complex;
	} else{
		std::cout<<"IncrementalPersistence does not support complexes\n";
		return nextPivots;
	}

	//Iterate over columns to reduce in reverse order
	for(auto columnIndexIter = edges.begin(); columnIndexIter != edges.end(); columnIndexIter++){
		simplexNodePointer simplex = (*columnIndexIter);	//The current simplex

		//Only need to test hash for equality
		//Not a pivot -> need to reduce
		if(it == pivots.end() || (*it)->hash != simplex->hash){

			//Get all cofacets using emergent pair optimization
			std::vector<nodeType *> faceList = (mode == "homology" ? complex->getAllFacets(simplex, saveVertices, dimension) \
			: (inData.complex->simplexType == "alphaComplex" ? inData.complex->getAllDelaunayCofacets_basePointer(simplex) \
			: complex->getAllCofacets(simplex, pivotPairs, true, saveVertices, dimension))); 

			std::vector<simplexNodePointer> columnV;	//Reduction column of matrix V
			columnV.push_back(simplex); //Initially V=I -> 1's along diagonal

			//Build a heap using the coface list to reduce and store in V
			std::make_heap(faceList.begin(), faceList.end(), compStruct);

			while(true){
				nodeType* pivot;

				while(!faceList.empty()){
					pivot = faceList.front();

					//Rotate the heap
					std::pop_heap(faceList.begin(), faceList.end(), compStruct);
					faceList.pop_back();

					if(!faceList.empty() && pivot->hash == faceList.front()->hash){ //Coface is in twice -> evaluates to 0 mod 2
						delete pivot;
						delete faceList.front();

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
				} else if(pivotPairs.find(pivot->hash) == pivotPairs.end()){ //Column cannot be reduced
					pivotPairs.insert({pivot->hash, simplex});
					nextPivots.push_back(std::shared_ptr<nodeType>(pivot));

					std::sort(columnV.begin(), columnV.end(), compStruct);
					auto it = columnV.begin();
					while(it != columnV.end()){
						if((it+1) != columnV.end() && (*it)->hash == (*(it+1))->hash) ++it;
						else v[simplex].push_back(*it);
						++it;
					}

					if(recordIntervals && simplex->weight != pivot->weight){
						//If we haven't saved the vertices in each simplex, we need to recover them now to determine the boundary points
						if(!saveVertices && inData.complex->simplexType != "alphaComplex"){
							for(simplexNodePointer simplex : v[simplex]){
								simplex->simplex = complex->getVertices(simplex->hash, dimension + 1, nPts);
							}
						}

						bettiBoundaryTableEntry des = { dimension, std::min(pivot->weight, simplex->weight), std::max(pivot->weight, simplex->weight), this->ut.extractBoundaryPoints(v[simplex]) };
						inData.bettiTable.push_back(des);
					}

					//Don't delete the first entry because that is converted to a smart pointer and stored as a pivot
					for(int i=1; i<faceList.size(); i++) delete faceList[i];

					break;
				} else{ 
					//Reduce the column of R by computing the appropriate columns of D by enumerating cofacets
					for(simplexNodePointer simp : v[pivotPairs[pivot->hash]]){
						columnV.push_back(simp);
						std::vector<nodeType *> faces = (mode == "homology" ? complex->getAllFacets(simp, saveVertices, dimension) : (inData.complex->simplexType == "alphaComplex" ? inData.complex->getAllDelaunayCofacets_basePointer(simp) : complex->getAllCofacets(simp, pivotPairs, false, saveVertices, dimension)));
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
//	IncrementalPersistence: For computing the persistence pairs from simplicial complex:
//		1. See Bauer-19 for algorithm/description
template<typename nodeType>
void incrementalPersistence<nodeType>::runPipe(pipePacket<nodeType> &inData){
	simplexArrayList<nodeType>* complex;	
	
	if(inData.complex->simplexType == "simplexArrayList" || inData.complex->simplexType == "alphaComplex" ){
		complex = (simplexArrayList<nodeType>*) inData.complex;
	} else{
		std::cout<<"IncrementalPersistence does not support complexes other than simplexArrayList and alphaComplex\n";
		return;
	}

	//Get the set of all points
	auto e = inData.complex->simplexType == "alphaComplex" ? complex->getdelaunayDimEdges(0) : complex->getDimEdges(0);
	//Convert the set to a vector
	std::vector<std::shared_ptr<nodeType>> edges = std::vector<std::shared_ptr<nodeType>>(e.begin(), e.end());
	//Initialize the binomial table
	complex->initBinom();
	//Get the next dimension (edges)
	edges = inData.complex->simplexType == "alphaComplex" ? complex->expanddelaunayDimension(1) : complex->expandDimension(edges);

	nPts = e.size();
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
	//TODO - POSSIBLY FIX
	std::vector<std::shared_ptr<nodeType>> pivots; //Store identified pivots
	unsigned mstSize = 0;
	unsigned nPts = inData.workData.size();

	unionFind uf(nPts);

	for(auto edgeIter = edges.begin(); edgeIter != edges.end(); edgeIter++){
		std::set<unsigned>::iterator it = (*edgeIter)->simplex.begin();

		//Find which connected component each vertex belongs to
		//	Use a hash map to track insertions for streaming or sparse indices
		unsigned v1 = *it;
		int c1 = uf.find(v1);

		it++;
		unsigned v2 = *it;
		int c2 = uf.find(v2);

		//Edge connects two different components -> add to the MST
		if(c1 != c2){
			uf.join(c1, c2);
			mstSize++;

			std::shared_ptr<nodeType> temp = std::make_shared<nodeType>(nodeType((*edgeIter)->simplex, (*edgeIter)->weight));
			temp->hash = v1 + v2*(v2-1)/2;
			pivots.push_back(temp);

			bettiBoundaryTableEntry des = { 0, 0, (*edgeIter)->weight, temp->simplex };
			inData.bettiTable.push_back(des);
		}

		//Check if we've filled our MST and can break
		if(mstSize >= edges.size()-1) break;
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

	for(unsigned d = 1; d < dim && !edges.empty(); d++){
		// //If d=1, we have already expanded the points into edges
		// //Otherwise, we need to generate the higher dimensional edges (equivalent to simplexList[d])
		if (d != 1)	edges = inData.complex->simplexType == "alphaComplex" ? complex->expanddelaunayDimension(d) : complex->expandDimension(edges, saveVertices, d);

		pivots = incrementalByDimension(inData, edges, pivots, d, sortReverseLexicographic(), "cohomology", !involuted);

		//To recover the representative cycles from the cocycles, we compute homology on just the pivot columns
		if(involuted){
			incrementalByDimension(inData, pivots, std::vector<std::shared_ptr<nodeType>>(), d, sortLexicographic(), "homology", true);
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
template<typename nodeType>
void incrementalPersistence<nodeType>::outputData(pipePacket<nodeType> &inData){
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
template<typename nodeType>
bool incrementalPersistence<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
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

	this->configured = true;
	this->ut.writeDebug("incrementalPersistence","Configured with parameters { dim: " + configMap["dimensions"] + ", complexType: " + configMap["complexType"] + ", eps: " + configMap["epsilon"]);
	this->ut.writeDebug("incrementalPersistence","\t\t\t\tdebug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

template class incrementalPersistence<simplexNode>;
template class incrementalPersistence<alphaNode>;
template class incrementalPersistence<witnessNode>;
