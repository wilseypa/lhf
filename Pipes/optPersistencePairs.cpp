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

std::pair<std::set<unsigned>,std::set<unsigned>> optPersistencePairs::getRankNull(std::vector<std::set<unsigned>> boundaryMatrix){
	//Perform column echelon reduction; basically the inverse of the RREF
	//	return column index of pivots (ranks)
	//	return row index of maximal faces for clearing
	
	utils ut;
		
	std::set<unsigned> retPivots;
	std::set<unsigned> retFaces;
	
	//Step through each column and search for pivot sequentially
	for(unsigned i = 0; i < boundaryMatrix.size(); i++){
		
		//Step through each simplex
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If the vector has j in the target set
			if(boundaryMatrix[j].size() > 0 && *(boundaryMatrix[j].begin()) == i){
				retPivots.insert(j);
	
				//Create a temp set for XOR
				std::set<unsigned> tempColumn = boundaryMatrix[j];
				
				//Record the lowest maximal face
				unsigned lastMaximal = *boundaryMatrix[j].begin();
				
				for(auto mf = boundaryMatrix[j].begin(); mf != boundaryMatrix[j].end(); mf++){
					lastMaximal = *mf;
				}
				if(retFaces.find(lastMaximal) == retFaces.end() && lastMaximal != *boundaryMatrix[j].begin())
					retFaces.insert(lastMaximal);
				
				
				//Iterate through remaining sets and XOR with our set
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j].size() > 0 && *(boundaryMatrix[j].begin()) == i){
						//XOR the remaining rows
						boundaryMatrix[j] = ut.setXOR(tempColumn, boundaryMatrix[j]);
					}
					j += 1;
				}				
			}
		}
	}
	std::cout << "RetPivots: " << retPivots.size() << "\t";
	ut.print1DVector(retPivots);
	std::cout << "RetFaces: " << retFaces.size() << "\t";
	ut.print1DVector(retFaces);
	
	return std::make_pair(retPivots, retFaces);
	
	
	
	
}

std::pair<std::set<unsigned>,std::set<unsigned>> optPersistencePairs::getRankNull(std::vector<std::vector<unsigned>> boundaryMatrix){
	
	//Perform column echelon reduction; basically the inverse of the RREF
	//	return column index of pivots (ranks)
	//	return row index of maximal faces for clearing
		
	std::set<unsigned> retPivots;
	std::set<unsigned> retFaces;
	
	//Step through each column and search for a one (pivot)
	for(unsigned i = 0; i < boundaryMatrix[0].size(); i++){
		
		//Step through each row
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If the vector has a 1 in the target column
			if(boundaryMatrix[j][i] == 1){
				retPivots.insert(j);
	
				//Create a temp matrix for XOR
				std::vector<unsigned> tempColumn = boundaryMatrix[j];
				
				unsigned lastMaximal = j;
				//Record the lowest maximal face
				for(unsigned mf = i + 1; mf < boundaryMatrix[0].size(); mf++){
					if(boundaryMatrix[j][mf] == 1){ //&& retFaces.find(mf) == retFaces.end() ){
						lastMaximal = mf;
					}
				}
				if(retFaces.find(lastMaximal) == retFaces.end())
					retFaces.insert(lastMaximal);
				
				
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
	std::cout << "RetPivots: " << retPivots.size() << "\t";
	ut.print1DVector(retPivots);
	std::cout << "RetFaces: " << retFaces.size() << "\t";
	ut.print1DVector(retFaces);
	return std::make_pair(retPivots, retFaces);
}

std::vector<std::set<unsigned>> optPersistencePairs::createBoundarySets(std::vector<std::vector<indSimplexTree::graphEntry>> ge, int d, std::set<unsigned> pivots, pipePacket p){
	
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();	
	
	std::vector<indSimplexTree::graphEntry> nChain;
	std::vector<indSimplexTree::graphEntry> pChain = ge[d];
	
	
	if(d == 0){
		std::vector<std::set<unsigned>> a;
		return a;
	} else
		nChain = ge[d-1];
		
	std::set<unsigned> a;
	std::vector<std::set<unsigned>> tempBoundary(pChain.size(), a);
	
	
	unsigned curPivot = -1;
	if(!pivots.empty()){
		curPivot = *(pivots.begin());
		//Create the columns (pChain), indexed by the rows in our pivot table
		for(unsigned i = 0; i < pChain.size(); i++){
			if(i == curPivot){
				if(!pivots.empty()){
					pivots.erase(pivots.begin());
					if(!pivots.empty())
						curPivot = *(pivots.begin());
				}				
			} else {
				tempBoundary[i] = pChain[i].getFaces(p.workData.complex);
			}
			
			
		}
			
		std::cout << std::endl;
	}else {	;
		for(unsigned i = 0; i < pChain.size(); i++){
			//Create the rows (nChain)
			tempBoundary[i] = pChain[i].getFaces(p.workData.complex);
		} 
		
	}
	return tempBoundary;
	
}


std::vector<std::vector<unsigned>> optPersistencePairs::createBoundaryMatrix(std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges, int d, std::set<unsigned> pivots){
	
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	
	//Setup p-chains and n-chains
	std::vector<std::pair<std::set<unsigned>,double>> nChain;
	std::vector<std::pair<std::set<unsigned>,double>> pChain = edges[d];

	if(d == 0){
		std::vector<std::vector<unsigned>> a;
		return a;
	}
	else
		nChain = edges[d-1];
	
	//Allocate the entire vector in a single step to reduce resizing during creation
	std::vector<std::vector<unsigned>> tempBoundary (pChain.size(), std::vector<unsigned>(nChain.size(), 0));
	
	unsigned curPivot = -1;
	if(!pivots.empty()){
		curPivot = *(pivots.begin());
		
		//Create the columns (pChain), indexed by the rows in our pivot table
		for(unsigned i = 0; i < pChain.size(); i++){
			if(i == curPivot){	
				
				if(!pivots.empty()){
					pivots.erase(pivots.begin());
					if(!pivots.empty())
						curPivot = *(pivots.begin());
				}
			} else {
						
				//Create the rows (nChain)
				for(unsigned j = 0; j < nChain.size(); j++){
					tempBoundary[i][j] = (checkFace(nChain[j].first, pChain[i].first));
				}
			}
		} 
	}else {	
		for(unsigned i = 0; i < pChain.size(); i++){
				
			//Create the rows (nChain)
			for(unsigned j = 0; j < nChain.size(); j++){
				tempBoundary[i][j] = (checkFace(nChain[j].first, pChain[i].first));
			}
		} 
		
	}
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();

	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;

	//Output the time and memory used for this pipeline segment
	std::cout << "Boundary Matrix (d=" << d << ") created in: " << (elapsed.count()/1000.0) << " seconds (physical time)" << std::endl;

	if(debug == 1){
		std::cout << std::endl << "_____BOUNDARY______" << std::endl;
		for(auto z : tempBoundary){
			std::cout << "  ";
			for(auto y : z)
				std::cout << y << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;
	}
		
		

	
	return tempBoundary;
}


// runPipe -> Run the configured functions of this pipeline segment
//
//	persistencePairs: For computing the persistence pairs from simplicial complex:
//		1. See Zomorodian-05 for algorithm/description of tArray
//		
pipePacket optPersistencePairs::runPipe(pipePacket inData){
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	unsigned pivotOffset = 0;
	std::set<unsigned> pivots;
	std::string bettis = "curIndex,index,dim,birth,death\n";
	auto edgeEndTime = std::chrono::high_resolution_clock::now();
	auto edgeStartTime = std::chrono::high_resolution_clock::now();
	
	std::set<unsigned> allPivots;
	std::set<unsigned> maxFaces;
	
	
	if(alterPipe){
		std::vector<std::vector<indSimplexTree::graphEntry>> indGraph = inData.workData.complex->getIndexEdges(maxEpsilon);
		edgeEndTime = std::chrono::high_resolution_clock::now();
		
		for(int d = dim; d > 0; d--){
			
			auto r = createBoundarySets(indGraph, d, maxFaces, inData);
			
			std::cout << "Getting Rank Nulls : " << r.size() << std::endl;
			//Compute the pivots of the current boundary matrix
			auto pivot_maxface_pair = getRankNull(r);
		
			pivots = pivot_maxface_pair.first;
			maxFaces = pivot_maxface_pair.second;
			
			pivotOffset -= inData.workData.complex->dimCounts[d];
			//Store the current pivots into allPivots (with offset) to compute tArray
			for(auto it = pivots.begin(); it != pivots.end() ; it++){
				unsigned val = *it;
				allPivots.insert(val + pivotOffset);
			}
			
			for(auto it = maxFaces.begin(); it != maxFaces.end() ; it++){
				unsigned val = *it;
			}
			
			
			//Iterate to next dimension; use the pivots to clear rows of the next dimension, use the boundary (later) 
			//		to compute the tArray
			
			
			
		}
		
		
	}
	
	else{
		std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges = inData.workData.complex->getAllEdges(maxEpsilon);
		edgeEndTime = std::chrono::high_resolution_clock::now();
		
		for(auto i : edges){
			pivotOffset += i.size();
		}
		
		
		std::vector<std::vector<unsigned>> temp;
		std::vector<std::vector<std::vector<unsigned>>> ret;
		for(int i = 0; i <= dim; i++){
			ret.push_back(temp);
		}
		

		//Iterate each dimension, only need coefficients (V)
		for(int d = dim; d > 0; d--){
			
			//First, construct the current boundary matrix (by dimension); clear pivots (if they exist)
			ret[d] = createBoundaryMatrix(edges, d, maxFaces);
				
			//Compute the pivots of the current boundary matrix
			auto pivot_maxface_pair = getRankNull(ret[d]);
		
			pivots = pivot_maxface_pair.first;
			maxFaces = pivot_maxface_pair.second;
			
			pivotOffset -= edges[d].size();
			//Store the current pivots into allPivots (with offset) to compute tArray
			for(auto it = pivots.begin(); it != pivots.end() ; it++){
				unsigned val = *it;
				allPivots.insert(val + pivotOffset);
			}
			
			for(auto it = maxFaces.begin(); it != maxFaces.end() ; it++){
				unsigned val = *it;
			}
			
			
			//Iterate to next dimension; use the pivots to clear rows of the next dimension, use the boundary (later) 
			//		to compute the tArray
		}	
				
				
		//Create the tArray from the identified pivots
		
		
		std::vector<std::pair<double,double>> temp2;
		std::vector<std::vector<std::pair<double,double>>> retarray;
		for(int i = 0; i < dim; i++){
			retarray.push_back(temp2);
		}
		
		
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
		if(!allPivots.empty()){
			curPivot = *(allPivots.begin());
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
					if(curDim == 2 && tArray[iter].marked && tArray[iter].death < 0 && tArray[curIndex].simplex.size() != kSimplices[iter].size())
						break;
					
					else if(tArray[iter].marked && tArray[iter].death < 0 \
						&& tArray[curIndex].simplex.size() != kSimplices[iter].size() \
						&& ut.setIntersect(kSimplices[iter], tArray[curIndex].simplex, false).size() == tArray[curIndex].simplex.size() - 1)
							break;
				}
				if(iter != 0){
					if(tArray[iter].death < 0 ){
						tArray[iter].death = kWeights[curIndex];
						
						if(tArray[iter].death != tArray[iter].birth)
							bettis += std::to_string(tArray[curIndex].simplex.size() - 2) + "," + std::to_string(tArray[iter].birth) + "," + std::to_string(kWeights[curIndex]) + "\n";
					}
				}
				
				if(!allPivots.empty()){
					allPivots.erase(allPivots.begin());
					if(!allPivots.empty()){
						curPivot = *(allPivots.begin());
					}
				}
			}
		}
	

		for(int t = 0; t < kSimplices.size(); t++){
			if(tArray[t].marked && tArray[t].death == -1 && tArray[t].simplex.size() <= dim){
				retarray[tArray[t].simplex.size() - 1].push_back(std::make_pair(kWeights[t], maxEpsilon));
				bettis += std::to_string(tArray[t].simplex.size() - 1) + "," + std::to_string(tArray[t].birth) + "," + std::to_string(maxEpsilon) + "\n";
			}
		}
	}
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();
	
	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
	std::chrono::duration<double, std::milli> edgeElapsed = edgeEndTime - edgeStartTime;
	
	//Output the time and memory used for this pipeline segment
	std::cout << "Bettis executed in " << (elapsed.count()/1000.0) << " seconds (physical time)" << std::endl;
	std::cout << "Edges executed in " << (edgeElapsed.count()/1000.0) << " seconds (physical time)" << std::endl;
	
	//Print the bettis	
	//std::cout << bettis << std::endl;
		
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
	
	file << "ti,Dim,Marked,Birth,Death,Simplex\n";
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
	
	pipe = configMap.find("complexType");
	if(pipe != configMap.end() && configMap["complexType"] == "indSimplexTree")
		alterPipe = true;
	
	return true;
}

