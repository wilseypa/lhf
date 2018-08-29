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
#include "utils.hpp"


std::vector<int> ranks;
std::vector<int> nRanks;
utils ut;

// basePipe constructor
bettiPipe::bettiPipe(){
	pipeType = "Betti";
	return;
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
	ut.print1DVector(face);
	ut.print1DVector(simplex);
	std::cout << "SIZE: " << ut.intersect(face,simplex,false).first.size() << std::endl;
	ut.print1DVector(ut.intersect(face,simplex,false).first);
	
	if(simplex.size() == 0)
		return 1;
	else if(ut.intersect(face,simplex,false).first.size() == face.size()){
		//ut.print1DVector(face);
		//ut.print1DVector(simplex);
		return 1;
	}
	else
		return 0;
}

// Reduce the binary boundary matrix
std::vector<std::vector<unsigned>> reduceBoundaryMatrix(std::vector<std::vector<unsigned>> boundaryMatrix){
	std::vector<std::vector<unsigned>> ret;
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
				//Along with a vector to clear out the row
				ret.push_back(boundaryMatrix[j]);
				auto tempRow = boundaryMatrix[j];
				auto clearRow = boundaryMatrix[j];
				clearRow[i] = 0;
				
				//Set up the temp row for both operations:
				for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
					tempRow[d] = tempRow[d] ^ clearRow[d];
					
				ut.print1DVector(tempRow);
				ut.print1DVector(clearRow);
				
				//Iterate through remaining vectors and XOR with our vector or the clearRow
				j += 1;
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j][i] == 1){
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] ^ tempRow[d];
					}
					else{
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] ^ clearRow[d];
					}
							
					j += 1;
				}			
				
				//Remove our vector from the boundaryMatrix and continue processing next column
				boundaryMatrix.erase(std::find(boundaryMatrix.begin(), boundaryMatrix.end(), tempRow));
			}
		}		
	}

	ranks.push_back(rank);
	nRanks.push_back(boundaryMatrix[0].size() - rank);
	std::cout << "COUNT:" << ret.size() << std::endl;
	return boundaryMatrix;
}

// Get the boundary matrix from the edges
std::vector<std::vector<unsigned>> boundaryMatrix(std::vector<std::vector<unsigned>> edges, int dim){
	std::vector<std::vector<unsigned>> ret;
	std::vector<std::vector<unsigned>> boundary;	
	
	//if(dim <= 1)
	//	return {{0}};
	
	//Create the boundary matrix from chains
	std::vector<std::vector<unsigned>> nChain = nSimplices(dim+2, edges);
	std::vector<std::vector<unsigned>> pChain = nSimplices(dim+1, edges);
	std::vector<unsigned> bDim;
	
	std::cout << dim << "-dimensional simplices: " << edges.size() << std::endl;
	
	std::cout << "pre: " << nChain.size() << " " << pChain.size() << std::endl;
	
	if (nChain.size() == 0)
		return {{0}};
	
	//
	if (pChain.size() == 0){
		for(unsigned j = 0; j < nChain.size(); j++){
			bDim.push_back(1);
		}
		boundary.push_back(bDim);
		bDim.clear();
	}
	else{
		for(unsigned i = 0; i < pChain.size(); i++){
			for(unsigned j = 0; j < nChain.size(); j++){
				bDim.push_back(checkFace(pChain[i], nChain[j]));
			}
			boundary.push_back(bDim);
			bDim.clear();
		}	
	}
		
	std::cout << "\tFinished boundary matrix -> " << boundary.size() << " x " <<boundary[1].size() << std::endl << std::endl;
	
	for(unsigned d = 0; d < boundary.size(); d++)
		ut.print1DVector(boundary[d]);
	std::cout << std::endl << std::endl;
	
	//Reduce the boundary matrix
	return reduceBoundaryMatrix(boundary);
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket bettiPipe::runPipe(pipePacket inData){
	std::vector<std::vector<unsigned>> local_edges;	
	std::vector<std::vector<std::vector<unsigned>>> allBoundaries;
	
	
	local_edges = inData.workData.complex->getEdges(0,0);
	
	std::cout << local_edges.size() << std::endl;
	std::cout << inData.workData.complex->simplexType << std::endl;
	
	for(int d = 0; d < dim; d++){		
		std::vector<std::vector<unsigned>> boundary = boundaryMatrix(local_edges, d);
		std::cout << "\tReduced boundary matrix -> " << boundary.size() ;
		std::cout << " x " <<boundary[1].size() << std::endl;
		
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
	if(pipe != configMap.end())
		dim = std::atoi(configMap["dimensions"].c_str());
	else return false;
	
	pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	else return false;
	
	return true;
}

