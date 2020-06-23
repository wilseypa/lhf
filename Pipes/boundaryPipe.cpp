/*
 * boundaryPipe hpp + cpp extend the basePipe class for calculating the 
 * minimum boundary from a distance matrix, used for upscaling the 
 * data for recomputation
 * 
 */

#include <string>
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
#include "boundaryPipe.hpp"

// basePipe constructor
boundaryPipe::boundaryPipe(){
	pipeType = "Boundary";
	return;
}

//Filter and return simplices of a specified dimension
std::vector<std::vector<unsigned>> boundaryPipe::nSimplices(double epsilon, unsigned n, std::vector<std::pair<double,std::vector<unsigned>>> complex){
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
int boundaryPipe::checkFace(std::vector<unsigned> face, std::vector<unsigned> simplex){
	
	if(simplex.size() == 0)
		return 1;
	else if(ut.symmetricDiff(face,simplex,false).size() == 1){
		return 1;
	}
	else
		return 0;
}

// Reduce the binary boundary matrix
std::pair<std::vector<std::vector<unsigned>>,std::pair<int,int>> boundaryPipe::reduceBoundaryMatrix(std::vector<std::vector<unsigned>> boundaryMatrix){
	std::vector<std::vector<unsigned>> ret;
	int rank = 0;
	
	if(boundaryMatrix.size() <= 0)
		return std::make_pair(ret,std::make_pair(0,0));
		
	//Step through each column and search for a 1 in that column
	for(unsigned i = 0; i < boundaryMatrix[0].size(); i++){
		//Step through each vector
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If the vector has a 1 in the target column
			if(boundaryMatrix[j][i] == 1){
				rank += 1;
				
				//Add the vector to our returned matrix
				auto tempRow = boundaryMatrix[j];
				ret.push_back(tempRow);
					
				//Remove our vector from the boundaryMatrix and continue processing next column
				boundaryMatrix.erase(boundaryMatrix.begin() + j);
				
				//Iterate through remaining vectors and XOR with our vector
				//j += 1;
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j][i] == 1){
						//XOR the remaining rows
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] ^ tempRow[d];
					}
					j += 1;
				}		
				for(auto a : ret){
					if (a[i] == 1){
						for(unsigned d = 0; d < ret.size(); d++)
							a[d] = a[d] ^ tempRow[d];
					}					
				}
			}
		}		
	}

	if(boundaryMatrix.size() > 0){
		for (auto a : boundaryMatrix)
			ret.push_back(a);
	}
	
	return std::make_pair(ret, std::make_pair(rank, (ret[0].size() - rank)));
}

// Get the boundary matrix from the edges
std::pair<std::vector<std::vector<unsigned>>,std::pair<int,int>> boundaryPipe::boundaryMatrix(std::vector<std::vector<unsigned>> nChain, std::vector<std::vector<unsigned>> pChain){
	std::vector<std::vector<unsigned>> ret;
	std::vector<std::vector<unsigned>> boundary;	
	
	//Create the boundary matrix from chains
	std::vector<unsigned> bDim;
	
	//std::cout << "pre: " << nChain.size() << " " << pChain.size() << std::endl;
	
	if (nChain.size() == 0)
		return std::make_pair(ret,std::make_pair(0,pChain.size()));
	
	//
	for(unsigned i = 0; i < nChain.size(); i++){
		std::vector<unsigned> bDim;
		
		//Create the columns (pChain)
		for(unsigned j = 0; j < pChain.size(); j++){
			bDim.push_back(checkFace(nChain[i], pChain[j]));
		}
		boundary.push_back(bDim);
	}	

	auto a = reduceBoundaryMatrix(boundary);
	
	
	
	return a;
}


std::vector<std::vector<unsigned>> boundaryPipe::extractBoundaries(std::vector<std::vector<unsigned>> edges, std::vector<std::vector<unsigned>> boundaryMatrix, int nullity){
	std::vector<std::vector<unsigned>> boundaries;
	
	//Iterate through each boundary 
	//	nullity = (dim - rank)
	//	Aggregate edges that form the constituent boundary
	
	// First, the columns (i to nullity), but they aren't sorted
	// Extract the nullity number of columns, these represent the loops
	//		Check each column for a leading 1
	for(int i = 0; i < boundaryMatrix[0].size(); i++){
		
		std::vector<unsigned> currentBoundary;		// Stores current boundary						// 
		
		//For each vector (row) in the boundary matrix 
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If there is a 1 in the row
			if(boundaryMatrix[j][i] == 1){
				
				for(auto edge : edges[j])
					if(std::find(currentBoundary.begin(), currentBoundary.end(), edge) == currentBoundary.end())
						currentBoundary.push_back(edge);
			}
		}
		
		std::cout << i << "\t" << nullity << std::endl;
		
		if(currentBoundary.size() > 1){
			std::cout << "Found Boundary\t";
			ut.print1DVector(currentBoundary);
			
			if(boundaries.size() == 0){
				boundaries.push_back(currentBoundary);
			} else {
				for(auto a : boundaries){
					if(ut.setIntersect(currentBoundary, a, false).size() == 0)
						boundaries.push_back(currentBoundary);
				}
			}
		}
		
	}
	std::cout << "ret" << std::endl;
	return boundaries;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//
pipePacket boundaryPipe::runPipe(pipePacket inData){
	
	if(dim < 2)
		return inData;
		
	int count = 0;
		
	std::vector<bettiBoundaryTableEntry> tempBetti;
	
	
	//Extract the higher-dimensional boundary points by using the map created from fastPersistence
	for(auto bet : inData.bettiTable){
		if(bet.bettiDim > 0){
			std::set<unsigned> totalBoundary;
			for(auto index : bet.boundaryPoints){
				totalBoundary.insert(index);
			}
			count++;
			tempBetti.push_back({bet.bettiDim, bet.birth, bet.death, totalBoundary, bet.boundary});
			
		} else {
			tempBetti.push_back(bet);
		}
	}
	
	std::cout << "Boundaries updated: " << count << std::endl;
	
	return inData;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void boundaryPipe::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_bettis_output.csv");
	
	file << inData.bettiOutput;
	
	file.close();
	return;
}

// configPipe -> configure the function settings of this pipeline segment
bool boundaryPipe::configPipe(std::map<std::string, std::string> configMap){
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
	
	configured = true;
	ut.writeDebug("boundary","Configured with parameters { dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	
	return true;
}

