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
#include "involutedPersistence.hpp"
#include "utils.hpp"

// basePipe constructor
involutedPersistence::involutedPersistence(){
	pipeType = "InvolutedPersistence";
	return;
}

void involutedPersistence::setupSimplices(std::vector<simplexNode*> simplices, unsigned d){
	edges = simplices;
	dim = d;
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//	involutedPersistence: For computing the persistence pairs from simplicial complex:
//		1. See Bauer-19 for algorithm/description
void involutedPersistence::runPipe(pipePacket &inData){
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();

	std::sort(edges.begin(), edges.end(), sortLexicographic());

	std::unordered_map<simplexNode*, std::vector<simplexNode*>> v;				//Store only the reduction matrix V and compute R implicity
	std::unordered_map<long long, simplexNode*> pivotPairs;						//For each pivot, which column has that pivot

	//Iterate over columns to reduce in reverse order
	for(auto columnIndexIter = edges.begin(); columnIndexIter != edges.end(); columnIndexIter++){
		simplexNode* simplex = (*columnIndexIter);		//The current simplex

		//Get all cofacets using emergent pair optimization
		std::vector<simplexNode*> faceList = inData.complex->getAllFacets(simplex);

		std::vector<simplexNode*> columnV;	//Reduction column of matrix V
		columnV.push_back(simplex); //Initially V=I -> 1's along diagonal

		//Build a heap using the coface list to reduce and store in V
		std::make_heap(faceList.begin(), faceList.end(), sortLexicographic());

		while(true){
			simplexNode* pivot;

			while(!faceList.empty()){
				pivot = faceList.front();

				//Rotate the heap
				std::pop_heap(faceList.begin(), faceList.end(), sortLexicographic());
				faceList.pop_back();

				if(!faceList.empty() && pivot->hash == faceList.front()->hash){ //Coface is in twice -> evaluates to 0 mod 2
					delete pivot;
					delete faceList.front();

					//Rotate the heap
					std::pop_heap(faceList.begin(), faceList.end(), sortLexicographic());
					faceList.pop_back();
				} else{

					faceList.push_back(pivot);
					std::push_heap(faceList.begin(), faceList.end(), sortLexicographic());
					break;
				}
			}

			if(faceList.empty()){ //Column completely reduced
				break;
			} else if(pivotPairs.find(pivot->hash) == pivotPairs.end()){ //Column cannot be reduced

				pivotPairs.insert({pivot->hash, simplex});

				std::sort(columnV.begin(), columnV.end());
				auto it = columnV.begin();
				while(it != columnV.end()){
					if((it+1) != columnV.end() && *it==*(it+1)) ++it;
					else v[simplex].push_back(*it);
					++it;
				}

				//Don't delete the first entry because that is converted to a smart pointer and stored as a pivot
				for(int i=1; i<faceList.size(); i++) delete faceList[i];

				if(simplex->weight != pivot->weight){
					bettiBoundaryTableEntry des = { dim, pivot->weight, simplex->weight, ut.extractBoundaryPoints(v[simplex]) };
					inData.bettiTable.push_back(des);
				}

				break;
			} else{ //Reduce the column of R by computing the appropriate columns of D by enumerating cofacets
				for(simplexNode* simp : v[pivotPairs[pivot->hash]]){
					columnV.push_back(simp);
					std::vector<simplexNode*> faces = inData.complex->getAllFacets(simp);
					faceList.insert(faceList.end(), faces.begin(), faces.end());
				}
				std::make_heap(faceList.begin(), faceList.end(), sortLexicographic());
			}
		}
	}

	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();

	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;

	//Output the time and memory used for this pipeline segment
	ut.writeDebug("persistence","Bettis executed in " + std::to_string(elapsed.count()/1000.0) + " seconds (physical time)");;

	return;
}



// outputData -> used for tracking each stage of the pipeline's data output without runtime
void involutedPersistence::outputData(pipePacket &inData){
	std::ofstream file;
	if(fnmod.size() > 0)
		file.open("output/"+pipeType+"_bettis_output"+fnmod+".csv");
	else
		file.open("output/" + pipeType + "_bettis_output.csv");

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
bool involutedPersistence::configPipe(std::map<std::string, std::string> &configMap){
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

	pipe = configMap.find("fn");
	if(pipe != configMap.end())
		fnmod = configMap["fn"];

	configured = true;
	ut.writeDebug("involutedPersistence","Configured with parameters { dim: " + configMap["dimensions"] + ", complexType: " + configMap["complexType"] + ", eps: " + configMap["epsilon"]);
	ut.writeDebug("involutedPersistence","\t\t\t\tdebug: " + strDebug + ", outputFile: " + outputFile + " }");

	return true;
}
