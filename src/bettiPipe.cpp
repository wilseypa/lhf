/*
 * bettiPipe hpp + cpp extend the basePipe class for calculating the 
 * betti numbers from a distance matrix
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include "bettiPipe.hpp"


std::vector<int> ranks;
std::vector<int> nRanks;

// basePipe constructor
bettiPipe::bettiPipe(){
	pipeType = "Betti";
	return;
}


// Find the intersect of two vectors
std::vector<unsigned> intersect(std::vector<unsigned> v1, std::vector<unsigned> v2){
	std::vector<unsigned> ret;
	
	//sort(v1.begin(), v1.end());
	//sort(v2.begin(), v2.end());
	
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(ret));
	
	return ret;
}

//Filter and return simplices of a specified dimension
std::vector<std::vector<unsigned>> nSimplices(unsigned n, std::vector<std::vector<unsigned>> complex){
	std::vector<std::vector<unsigned>> ret;
	
	for(auto v : complex){
		if(v.size() == n)
			ret.push_back(v);
	}
	
	return ret;
}

// Check if a face is a subset of a simplex
int checkFace(std::vector<unsigned> face, std::vector<unsigned> simplex){
	if(simplex.size() == 0)
		return 1;
	else if(intersect(face,simplex).size() == face.size())
		return 1;
	else
		return 0;
}

// Reduce the binary boundary matrix
std::vector<std::vector<int>> reduceBoundaryMatrix(std::vector<std::vector<int>> boundaryMatrix){
	std::vector<std::vector<int>> ret;
	int rank = 0;
	if(boundaryMatrix.size() <= 0)
		return boundaryMatrix;
		
	//Step through each column and search for a 1 in that column
	for(unsigned i = 0; i < boundaryMatrix[0].size(); i++){
		
		//Step through each vector
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If the vector has a 1 in the target column
			if(boundaryMatrix[j][i] == 1){
				rank += 1;
				
				//Add the vector to our returned matrix
				ret.push_back(boundaryMatrix[j]);
				auto tempRow = boundaryMatrix[j];
				
				//Iterate through remaining vectors and XOR with our vector
				j += 1;
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j][i] == 1)
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] ^ tempRow[d];
					j += 1;
				}			
				
				//Remove our vector from the boundaryMatrix and continue processing next column
				boundaryMatrix.erase(std::find(boundaryMatrix.begin(), boundaryMatrix.end(), tempRow));
			}
		}		
	}

	ranks.push_back(rank);
	nRanks.push_back(boundaryMatrix[0].size() - rank);
	return ret;
}

// Get the boundary matrix from the edges
std::vector<std::vector<int>> boundaryMatrix(std::vector<std::vector<unsigned>> edges, int dim){
	std::vector<std::vector<int>> ret;
	std::vector<std::vector<int>> boundary;	
	
	//if(dim <= 1)
	//	return {{0}};
	
	//Create the boundary matrix from chains
	std::vector<std::vector<unsigned>> nChain = nSimplices(dim, edges);
	std::vector<std::vector<unsigned>> pChain = nSimplices(dim-1, edges);
	std::vector<int> bDim;
	
	std::cout << "pre: " << nChain.size() << " " << pChain.size() << std::endl;
	
	if (nChain.size() == 0)
		return {{0}};
	
	//
	for(unsigned i = 0; i < nChain.size(); i++){
		if(pChain.size() == 0)
			bDim.push_back(1);
		for(unsigned j = 0; j < pChain.size(); j++){
			bDim.push_back(checkFace(pChain[j], nChain[i]));
		}
		boundary.push_back(bDim);
		bDim.clear();
	}	
		
	std::cout << "\tFinished boundary matrix -> " << boundary.size() << " x " <<boundary[1].size() << std::endl;
	
	//Reduce the boundary matrix
	return reduceBoundaryMatrix(boundary);
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket bettiPipe::runPipe(pipePacket inData){
	
	std::vector<std::vector<std::vector<int>>> allBoundaries;
	
	for(int d = 0; d < dim; d++){		
		std::vector<std::vector<int>> boundary = boundaryMatrix(inData.workData.edges, d);
		
		std::cout << "\tReduced boundary matrix -> " << boundary.size() << " x " <<boundary[1].size() << std::endl;
		
		std::cout << "\tDim: " << d << "\tRANKS: ";
		for(auto z : ranks)
			std::cout << z << " ";
		std::cout << std::endl;
	
		allBoundaries.push_back(boundary);
	}
	
	
	
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool bettiPipe::configPipe(std::map<std::string, std::string> configMap){
	
	
	auto pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	return true;
}

