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
	
	
	//Build buffers for each thread
	int n = threads;
	std::unordered_map<simplexNode*, std::vector<simplexNode*>> boundary;	//Store the boundary matrix
	//std::unordered_map<simplexNode*, simplexNode*> pivotPairs;				//For each pivot, which column has that pivot
	std::vector<std::pair<simplexNode*, std::vector<simplexNode*>>> columnsToReduce; //Columns in the jth range which need to be reduced
	
	simplexNode* first[n]; 	//First simplex in the ith range
	
	std::vector<unsigned> blocks;
	int nSimplices = 0; //Total number of simplices
	for(unsigned d = 0; d <= dim; d++){
		blocks.push_back(nSimplices);
		
		nSimplices += edges[d].size();
	}
	blocks.push_back(nSimplices);
	std::cout << "Total Simplices: " << nSimplices << std::endl;
	
	std::cout << "Blocks: \t";
	for(auto a : blocks)
		std::cout << a << "\t";
	std::cout << std::endl;
	
	std::cout << "d0: " << edges[0].size() << std::endl;
	std::cout << "d1: " << edges[1].size() << std::endl;
	std::cout << "d2: " << edges[2].size() << std::endl; 
	
	std::cout << "starting threads: " << threads << std::endl;
	
	std::vector<int> iter(threads, 0);
	std::vector<bool> complete(50, 0);
	
	//Process within each dimension
	for(unsigned d = 0; d <= dim; d++){
	
		//Within the dimension, the algorithm needs to process edges[d] columns; this defines our omp FOR clause
		//		The bounds of these points (i_min, i_max) - 
		//			i_min = blocks[d] ; i_max = blocks[d+1];
		//		We only consider these in the forms of columns of reduction to be processed with their cofaces
	
		int i_min = blocks[d];
		int i_max = blocks[d+1];
		
		std::cout << "Using bounds: " << i_min << "\t" << i_max << std::endl;
		
		//Reduce only to our necessary columns:
		
		columnsToReduce = inData.complex->simplexList[d];
		
		std::cout << "reducing columns: " << columnsToReduce.size() << std::endl;
		
		inData.complex->prepareCofacets(d);
		std::sort(pivots.begin(), pivots.end(), cmpByWeightDec());
		std::vector<simplexNode*>::iterator it = pivots.begin();

		std::vector<simplexNode*> nextPivots;	 					//Pivots for the next dimension
		std::unordered_map<simplexNode*, std::vector<simplexNode*>> v;				//Store only the reduction matrix V and compute R implicity
		std::unordered_map<simplexNode*, simplexNode*> pivotPairs;	//For each pivot, which column has that pivot
		
		unsigned curCoIndex = 0;
	
		#pragma omp parallel num_threads(threads)
		{
			int np = omp_get_thread_num();
			//
			//	For each column to reduce, which contains a list of cofacets and an index in the weighted list
			//		if the index before us (n-1) has been processed in the boolean vector, our reduced lowest cofacet is a pivot
			//			While the index before us hasn't been processed, we can reduce in parallel against the reduction matrix
			//				To do this, wait for the lowest cofacet to be a pivot, then reduce, and wait for the next lowest cofacet to be a pivot
			//
			//			
			//		Some Notes:
			//			-This is a np-independent algorithm; don't worry about n (# of threads) or np (current thread #) in this context
			//			
			//		Key Storage:
			//			-curCoIndex - stores the current index that is ready to be emitted (keeps sequential ordering of the columns)
			//			

			
			#pragma omp for schedule(dynamic) ordered
			for(unsigned coIndex = 0; coIndex < columnsToReduce.size(); coIndex++){ 
				
				auto simp = columnsToReduce[coIndex];
				simplexNode* simplex = simp.first;
				std::vector<simplexNode*> cofaceList = simp.second;
				
				if(pivotPairs[i].find(simplex) != pivotPairs[i].end()) continue;

				//Build a heap using the coface list to reduce and store in V
				if(cofaceList.empty()){
					cofaceList = inData.complex->getAllCofacets(simplex->simplex);
					std::make_heap(cofaceList.begin(), cofaceList.end(), cmpByWeightDec());
				}				
				
				
				//Loop through reducing the columns until we can push our pivot or clear
				while(curCoIndex != coIndex){
					
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
							(unreducedColumns).push_back(make_pair(simplex, cofaceList));
							break;
						} else if(pivotPairs.find(pivot) == pivotPairs.end()){ //Column cannot be reduced
							pivotPairs.insert({pivot, simplex});

							boundary[simplex] = cofaceList;

							//Needs to be synchronized or inserted into array
							if(simplex->weight != pivot->weight){
								bettiBoundaryTableEntry des = { simplex->simplex.size()-1, simplex->weight, pivot->weight, {}, cofaceList };
								inData.bettiTable.push_back(des);
							}

							break;
						} else{ //Reduce the column of R by computing the appropriate columns of D by enumerating cofacets

							auto cofaces = boundary[pivotPairs[pivot]];
							cofaceList.insert(cofaceList.end(), cofaces.begin(), cofaces.end());
							std::make_heap(cofaceList.begin(), cofaceList.end(), cmpByWeightDec());
						}
					}
				}			
				
				

						
			}
			
			
		}
		//END PARALLEL EXECUTION	
		
		std::cout << std::endl;
	}

	for(int i = 0; i < iter.size(); i++){
		std::cout << "Thread #" << i << " completed " << iter[i] << std::endl;
	}

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
