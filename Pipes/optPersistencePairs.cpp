/*
 * optPersistencePairs hpp + cpp extend the basePipe class for calculating the 
 * optimized persistence pairs numbers from a complex
 * 
 */

#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include "optPersistencePairs.hpp"


// basePipe constructor
optPersistencePairs::optPersistencePairs(){
	pipeType = "PersistencePairs";
	return;
}

//Filter and return simplices of a specified dimension
std::vector<std::vector<unsigned>> optPersistencePairs::nSimplices(double epsilon, unsigned n, std::vector<std::pair<double,std::vector<unsigned>>> complex){
	std::vector<std::vector<unsigned>> ret;
	for(auto v : complex){
		if(v.second.size() == n){
			if(v.first <= epsilon)
				ret.push_back(v.second);
		}
	}
	
	return ret;
}

// Check if a face is a subset of a simplex
int optPersistencePairs::checkFace(std::vector<unsigned> face, std::vector<unsigned> simplex){
	
	if(simplex.size() == 0)
		return 1;
	else if(ut.symmetricDiff(face,simplex,false).size() == 1){
		return 1;
	}
	else
		return 0;
}

// Check if a face is a subset of a simplex
int optPersistencePairs::checkFace(std::set<unsigned> face, std::set<unsigned> simplex){
	
	if(simplex.size() == 0)
		return 1;
	else if(ut.symmetricDiff(face,simplex,false).size() == 1){
		return 1;
	}
	else
		return 0;
}

std::set<unsigned> optPersistencePairs::getRankNull(std::vector<std::vector<unsigned>> boundaryMatrix, unsigned offset){
	
	//Perform column echelon reduction; basically the inverse of the RREF
	//	return offset index of pivots (ranks)
	
	std::set<unsigned> ret;
	
	//Step through each column and search for a one
	for(unsigned i = 0; i < boundaryMatrix[0].size(); i++){
		
		//Step through each row
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If the vector has a 1 in the target column
			if(boundaryMatrix[j][i] == 1){
				ret.insert(j);
	
				//Create a temp matrix for XOR
				std::vector<unsigned> tempColumn = boundaryMatrix[j];
				
				//Zero out our column and continue processing next column to clear any ones
				//for(int z = 0; z < boundaryMatrix.size(); z++){
				//	boundaryMatrix[j][i] = 0;
			//	}
				
				//Iterate through remaining vectors and XOR with our vector or the clearRow
				//j += 1;
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j][i] == 1){
						//XOR the remaining rows
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] ^ tempColumn[d];
					}
					j += 1;
				}				
			}
		}
	}
	
	return ret;
}

std::set<unsigned> optPersistencePairs::createBoundaryMatrix(std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges, int d){
	std::set<unsigned> pivots;
	std::set<unsigned> lastPivots;
	
	//Iterate through each dimensional graph...
	
	//Setup p-chains and n-chains
	std::vector<std::pair<std::set<unsigned>,double>> nChain;
	std::vector<std::pair<std::set<unsigned>,double>> pChain = edges[d];
	
	if(d == 0)
		nChain = {};
	else
		nChain = edges[d-1];
		
	//Allocate the entire vector in a single step to reduce resizing during creation
	std::vector<std::vector<unsigned>> tempBoundary (pChain.size(), std::vector<unsigned>(nChain.size(), 0));
		
	//Create the columns (pChain), indexed by the rows in our pivot table
	for(unsigned i = 0; i < pChain.size(); i++){
		//Create the columns (pChain)
		for(unsigned j = 0; j < nChain.size(); j++){
			tempBoundary[i][j] = (checkFace(nChain[j].first, pChain[i].first));
		}
	}
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	//Compute the pivots of the
	lastPivots = getRankNull(tempBoundary, 0);	
	
	for(auto a : lastPivots)
		pivots.insert(a);
		
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();

	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;

	//Output the time and memory used for this pipeline segment
	std::cout << "GetRankNull calculated for dim " << d << " executed in " << (elapsed.count()/1000.0) << " seconds (physical time)" << std::endl;

	std::cout << "\tPivots (" << pivots.size() << "): ";
	for (auto a : pivots)
		std::cout << a << ",";
	std::cout << std::endl;
	
	return pivots;
}


// runPipe -> Run the configured functions of this pipeline segment
//
//	persistencePairs: For computing the persistence pairs from simplicial complex:
//		1. See Zomorodian-05 for algorithm/description of tArray
//		
pipePacket optPersistencePairs::runPipe(pipePacket inData){
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges = inData.workData.complex->getAllEdges(maxEpsilon);
	std::string bettis = "curIndex,index,dim,birth,death\n";
	
	std::vector<std::pair<double,double>> temp;
	std::vector<std::vector<std::pair<double,double>>> ret;
	for(int i = 0; i < dim; i++){
		ret.push_back(temp);
	}
	
	//In basic terms
	//		For each dim, 0 <= dim <= k
	//			Iterate through total-ordered simplices, (by dim, weight)
	//				Call RemovePivotRows and retrieve simplex or null
	//					if null, mark the simplex (in T[i]) and continue
	//					if there is a ret, find that simplex in the T-array (T[i])
	//						Generate a persistence interval (weight of T[i], weight of T[j])
	//						Store j and d in T[i];
	//			Afterwards, find all marked rows, j, with no entry for T[j]
	//				Generate persistence interval for these, (weight of T[j], inf)
	//		Success!
	
	//Iterate each dimension, only need ranks
	for(int d = dim; d >= 0; d--){
		
		//Get pivots of current boundary matrix (d)
		std::set<unsigned> kPivots = createBoundaryMatrix(edges, d);
	
		
		//Flatten the edges into a single array
		std::vector<std::set<unsigned>> kSimplices;
		std::vector<double> kWeights;
		for(auto a : edges){		
			for(auto z : a){
				kSimplices.push_back(z.first);
				kWeights.push_back(z.second);
			}
		}
		
		unsigned curPivot = -1;
		if(!kPivots.empty()){
			curPivot = *(kPivots.begin());
		}
		
		
		
		for(int curIndex = 0; curIndex < kSimplices.size(); curIndex++){
			
			int curDim = kSimplices[curIndex].size();
			tArray.push_back(tArrayEntry_t());
			tArray[curIndex].simplex = kSimplices[curIndex];
			
			if(curDim == 1){
				tArray[curIndex].birth = kWeights[curIndex];
				tArray[curIndex].marked = true;
			}
			else if(curIndex != curPivot){
				tArray[curIndex].birth = kWeights[curIndex];
				tArray[curIndex].marked = true;
				
			} else {
				int iter = tArray.size();
				for(; iter >= 0; iter--){
					if(tArray[iter].marked && tArray[iter].death < 0 && tArray[curIndex].simplex.size() != kSimplices[iter].size() && ut.setIntersect(kSimplices[iter], tArray[curIndex].simplex, false).size() == tArray[curIndex].simplex.size() - 1)
						break;
				}
				if(iter != 0){
					if(tArray[iter].death < 0 ){
						tArray[iter].death = kWeights[curIndex];
						
						if(tArray[iter].death != tArray[iter].birth)
							bettis += std::to_string(tArray[curIndex].simplex.size() - 2) + "," + std::to_string(tArray[iter].birth) + "," + std::to_string(kWeights[curIndex]) + "\n";
					}
				}
				
				if(!kPivots.empty()){
					kPivots.erase(kPivots.begin());
					curPivot = *(kPivots.begin());
				}
			}
		}
	
	
		for(int t = 0; t < kSimplices.size(); t++){
			if(tArray[t].marked && tArray[t].death == -1 && tArray[t].simplex.size() < dim){
				ret[tArray[t].simplex.size()].push_back(std::make_pair(kWeights[t], maxEpsilon));
				bettis += std::to_string(0) + "," + std::to_string(tArray[t].birth) + "," + std::to_string(maxEpsilon) + "\n";
			}
		}
	}
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();
	
	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
	
	//Output the time and memory used for this pipeline segment
	std::cout << "Bettis executed in " << (elapsed.count()/1000.0) << " seconds (physical time)" << std::endl;
	
	//Print the bettis	
	std::cout << bettis << std::endl;
		
	inData.bettiOutput = bettis;
		
	return inData;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void optPersistencePairs::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_bettis_output.csv");
	
	file << inData.bettiOutput;
	
	file.close();
	
	file.open("output/tArray.csv");
	
	file << "ti\tDim\tMarked\tBirth\tDeath\tSimplex\n";
	for(auto tStruct : tArray){
		file << tStruct.ti << "," << tStruct.simplex.size()-1 << "," << tStruct.marked << "," << tStruct.birth << "," << tStruct.death << ",";
		for(auto index : tStruct.simplex)
			file << index << " ";
		file << "\n";
	}
	file.close();
	
	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool optPersistencePairs::configPipe(std::map<std::string, std::string> configMap){
	auto pipe = configMap.find("dimensions");
	if(pipe != configMap.end())
		dim = std::atoi(configMap["dimensions"].c_str());
	else return false;
	
	pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	else return false;
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return false;	
	
	pipe = configMap.find("twist");
	if(pipe != configMap.end())
		twist = configMap["twist"];
	else return false;
	
	return true;
}

