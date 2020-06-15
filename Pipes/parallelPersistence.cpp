/*
 * persistencePairs hpp + cpp extend the basePipe class for calculating the 
 * persistence pairs numbers from a complex
 * 
 */

#include "omp.h"
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
#include "parallelPersistence.hpp"
#include "utils.hpp"

// basePipe constructor
parallelPersistence::parallelPersistence(){
	pipeType = "ParallelPersistence";
	return;
}

struct cmpSimplices{
	bool operator()(simplexNode* a, simplexNode* b){
		if(a->simplex.size() == b->simplex.size()) return cmpByWeightDec()(a, b);
		return a->simplex.size() < b->simplex.size();
	}
};

// runPipe -> Run the configured functions of this pipeline segment
//
//	parallelPersistence: For computing the persistence pairs from simplicial complex:
//		1. See Bauer-19 for algorithm/description
pipePacket parallelPersistence::runPipe(pipePacket inData){
	
	//Get all edges for the simplexArrayList or simplexTree
	std::vector<std::set<simplexNode*, cmpByWeight>> edges = inData.complex->getAllEdges();

	if(edges.size() <= 1) return inData;
	
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel num_threads(threads)
	{
		int n = omp_get_num_threads();
		printf("Thread %d started\n", omp_get_thread_num());

		std::unordered_map<simplexNode*, std::vector<simplexNode*>> boundary[n];	//Store the boundary matrix
		std::unordered_map<simplexNode*, simplexNode*> pivotPairs[n];				//For each pivot, which column has that pivot
		std::vector<std::pair<simplexNode*, std::vector<simplexNode*>>> columnsToReduce[n]; //Columns in the jth range which need to be reduced
		simplexNode* first[n]; 	//First simplex in the ith range

		int nSimplices = 0; //Total number of simplices
		for(unsigned d = 0; d <= dim; d++) nSimplices += edges[d].size();

		int blockSize = nSimplices/n; //Create approximately equal size blocks
		std::vector<unsigned> blocks;
		for(int i = 0; i < n; i++) blocks.push_back(i*blockSize);
		blocks.push_back(nSimplices);

		int block = 0;
		unsigned i = 0;

		first[0] = *edges[0].rbegin();
		for(unsigned d = 0; d <= dim; d++){
			for(auto it = edges[d].rbegin(); it != edges[d].rend(); it++){
				if(i == blocks[block+1]){
					block++;
					first[block] = *it;
				}

				//Iterate over all the edges and assign each column to the correct thread
				columnsToReduce[block].push_back(make_pair(*it, std::vector<simplexNode*>()));
				i++;
			}
		}

		//Row range i (node i) and column range j
		//Iterate from (n-1, 0)
		// (n-2, 0), (n-1, 1)
		// (n-3, 0), (n-2, 1), (n-1, 2)
		// etc.
		for(int diff = n-1; diff >= 0; diff--){ //Difference between i and j
			for(int j = 0; j < n-diff; j++){
				int i = j + diff;

				//Columns that won't be reduced on this node -> must be sent to the next node
				std::vector<std::pair<simplexNode*, std::vector<simplexNode*>>> unreducedColumns;

				for(auto simp : columnsToReduce[j]){
					simplexNode* simplex = simp.first;
					std::vector<simplexNode*> cofaceList = simp.second;

					if(pivotPairs[i].find(simplex) != pivotPairs[i].end()) continue;

					//Build a heap using the coface list to reduce and store in V
					if(cofaceList.empty()){
						cofaceList = inData.complex->getAllCofacets(simplex->simplex);
						std::make_heap(cofaceList.begin(), cofaceList.end(), cmpByWeightDec());
					}

					while(true){
						simplexNode* pivot;
						while(!cofaceList.empty()){
							pivot = cofaceList.front();

							//Rotate the heap
							std::pop_heap(cofaceList.begin(), cofaceList.end(), cmpByWeightDec());
							cofaceList.pop_back();

							if(!cofaceList.empty() && pivot == cofaceList.front()){ //Coface is in twice -> evaluates to 0 mod 2
								
								//Rotate the heap
								std::pop_heap(cofaceList.begin(), cofaceList.end(), cmpByWeightDec());
								cofaceList.pop_back();
							} else{

								cofaceList.push_back(pivot);
								std::push_heap(cofaceList.begin(), cofaceList.end(), cmpByWeightDec());
								break;
							}
						}
						
						if(cofaceList.empty()){ //Column completely reduced
							break;
						} else if(cmpSimplices()(pivot, first[i])){ //Pivot is not in range i
							unreducedColumns.push_back(make_pair(simplex, cofaceList));
							break;
						} else if(pivotPairs[i].find(pivot) == pivotPairs[i].end()){ //Column cannot be reduced
							pivotPairs[i].insert({pivot, simplex});

							boundary[i][simplex] = cofaceList;

							if(simplex->weight != pivot->weight){
								bettiBoundaryTableEntry des = { simplex->simplex.size()-1, simplex->weight, pivot->weight, {}, cofaceList };
								inData.bettiTable.push_back(des);
							}

							break;
						} else{ //Reduce the column of R by computing the appropriate columns of D by enumerating cofacets

							auto cofaces = boundary[i][pivotPairs[i][pivot]];
							cofaceList.insert(cofaceList.end(), cofaces.begin(), cofaces.end());
							std::make_heap(cofaceList.begin(), cofaceList.end(), cmpByWeightDec());
						}
					}
				}
				//Send unreduced columns to the next node
				columnsToReduce[j] = unreducedColumns;
			}
		}
		
	}
	//END PARALLEL EXECUTION

	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();

	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;

	//Output the time and memory used for this pipeline segment
	ut.writeDebug("parallelPersistence","Bettis executed in " + std::to_string(elapsed.count()/1000.0) + " seconds (physical time)");;
		
	return inData;
}



// outputData -> used for tracking each stage of the pipeline's data output without runtime
void parallelPersistence::outputData(pipePacket inData){
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
bool parallelPersistence::configPipe(std::map<std::string, std::string> configMap){
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
	
	pipe = configMap.find("threads");
	if(pipe != configMap.end())
		threads = std::atoi(configMap["threads"].c_str());
	else return false;
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return false;	
	
	pipe = configMap.find("fn");
	if(pipe != configMap.end())
		fnmod = configMap["fn"];
		
	configured = true;
	ut.writeDebug("parallelPersistence","Configured with parameters { dim: " + configMap["dimensions"] + ", complexType: " + configMap["complexType"] + ", eps: " + configMap["epsilon"]);
	ut.writeDebug("parallelPersistence","\t\t\t\tdebug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}
