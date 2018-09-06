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
std::vector<std::vector<unsigned>> nSimplices(double epsilon, unsigned n, std::vector<std::pair<double,std::vector<unsigned>>> complex){
	std::vector<std::vector<unsigned>> ret;
	
	for(auto v : complex){
		
		//for(auto i : v.second)
		//	std::cout << i << " ";
		//std::cout << "--> " << v.second.size() << " " << v.first << " " << epsilon << std::endl << std::endl;
		
		
		if(v.second.size() == n){
			if(v.first <= epsilon)
				ret.push_back(v.second);
		}
	}
	
	//std::cout << "nSimplices : " << n << std::endl;
	//for(auto i : ret){
	//	for(auto p : i)
	//		std::cout << p << " ";
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	
	return ret;
}

// Check if a face is a subset of a simplex
int checkFace(std::vector<unsigned> face, std::vector<unsigned> simplex){
	//ut.print1DVector(face);
	//ut.print1DVector(simplex);
	//std::cout << "SIZE: " << ut.intersect(face,simplex,false).first.size() << std::endl;
	//ut.print1DVector(ut.intersect(face,simplex,false).first);
	
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
std::pair<std::vector<std::vector<unsigned>>,int> reduceBoundaryMatrix(std::vector<std::vector<unsigned>> boundaryMatrix){
	std::vector<std::vector<unsigned>> ret;
	int rank = 0;
	if(boundaryMatrix.size() <= 0)
		return std::make_pair(boundaryMatrix, 0);
		
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
				
				for(unsigned d = i+1; d < boundaryMatrix[0].size(); d++)
					clearRow[d] = tempRow[d];
					
				//Remove our vector from the boundaryMatrix and continue processing next column
				boundaryMatrix.erase(boundaryMatrix.begin() + j);
				
				//Iterate through remaining vectors and XOR with our vector or the clearRow
				//j += 1;
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j][i] == 1){
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] != tempRow[d];
					}
					//else{
						//for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							//boundaryMatrix[j][d] = boundaryMatrix[j][d] != clearRow[d];
					//}
							
					j += 1;
				}			
				
				for(auto a : ret){
					if (a[i] == 1){
						for(unsigned d = 0; d < ret.size(); d++)
							a[d] = a[d] != tempRow[d];
					}					
				}
			}
		}		
	}

	if(boundaryMatrix.size() > 0){
		for (auto a : boundaryMatrix)
			ret.push_back(a);
	}
	
	std::cout << "Rank: " << rank << std::endl;
	std::cout << "Nullity: " << (ret[0].size() - rank) << std::endl;
	ranks.push_back(rank);
	nRanks.push_back(ret.size() - rank);
	return std::make_pair(ret, rank);
}

// Get the boundary matrix from the edges
std::pair<std::vector<std::vector<unsigned>>,int> boundaryMatrix(std::vector<std::vector<unsigned>> nChain, std::vector<std::vector<unsigned>> pChain){
	std::vector<std::vector<unsigned>> ret;
	std::vector<std::vector<unsigned>> boundary;	
	
	
	//Create the boundary matrix from chains
	std::vector<unsigned> bDim;
	
	//std::cout << "pre: " << nChain.size() << " " << pChain.size() << std::endl;
	
	if (nChain.size() == 0)
		return std::make_pair(ret,0);
	
	//
	if (pChain.size() == 0){
		for(unsigned j = 0; j < nChain.size(); j++){
			bDim.push_back(1);
		}
		boundary.push_back(bDim);
		bDim.clear();
	}
	else{
		std::cout << nChain.size() << "\t" << pChain.size() << std::endl << std::endl;
		
		std::cout << "_______nChain______" << std::endl;
		for(auto i : nChain){
			for(auto p : i)
				std::cout << p << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;
		
		std::cout << "_______pChain______" << std::endl;
		for(auto i : pChain){
			for(auto p : i)
				std::cout << p << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
		
		
		//Create the rows (nChain)
		for(unsigned i = 0; i < pChain.size(); i++){
			//Create the columns (pChain)
			for(unsigned j = 0; j < nChain.size(); j++){
				bDim.push_back(checkFace(pChain[i], nChain[j]));
			}
			boundary.push_back(bDim);
			bDim.clear();
		}	
	}
		
	std::cout << "\tFinished boundary matrix -> " << boundary.size() << " x " <<boundary[1].size() << std::endl << std::endl;
	
	
	std::cout << std::endl << "ORIGINAL BOUNDARY" << std::endl;
	for(unsigned d = 0; d < boundary.size(); d++)
		ut.print1DVector(boundary[d]);
	std::cout << std::endl << std::endl;
	
	//Reduce the boundary matrix
	return reduceBoundaryMatrix(boundary);
}


std::vector<std::vector<unsigned>> extractBoundaries(std::vector<std::pair<double, std::vector<unsigned>>> edges, std::vector<std::vector<unsigned>> boundaryMatrix, int nullity){
	std::vector<std::vector<unsigned>> boundaries;
	
	//Iterate through each boundary 
	//	nullity = (dim - rank)
	//	Aggregate edges that form the constituent boundary
	for(int i = 0; i < nullity; i++){
		
		//store the current boundary
		std::vector<unsigned> currentBoundary;
		int j = 0;
		
		//
		for(auto vect : boundaryMatrix){
			
			if(vect[vect.size() - nullity] == 1){
				for(auto edge : edges[j].second){
					std::cout << edge << " ";
					currentBoundary.push_back(edge);
				}
			}
			j++;
			
		}
		
		for(auto edge : edges[edges.size() - nullity + i].second){
			std::cout << edge << " ";
			currentBoundary.push_back(edge);
		}
		std::cout << std::endl << std::endl;
		
		for (auto e : currentBoundary)
			std::cout << e << ",";
		std::cout << std::endl;
		
		boundaries.push_back(currentBoundary);
		
	}
	
	return boundaries;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//	bettiPipe: generate boundary matrix, determine constituent boundaries
//		edges that make up the boundary
//
//		For each dimension, 
//			Generate boundary matrix at each weight
//				Then, reduce the boundary matrix to RREF
//				Extract boundaries from RREF
//			Generate barcodes (lifespans) of each component
//				WITH boundaries that form the barcode
//		
//
pipePacket bettiPipe::runPipe(pipePacket inData){
	std::vector<std::vector<std::vector<unsigned>>> allBoundaries;

	
	auto local_edges = inData.workData.complex->getEdges(0,0);
	
	//std::cout << local_edges.size() << std::endl;
	//std::cout << inData.workData.complex->simplexType << std::endl;
	
	
	//Iterate through each dimension to build boundary matrix
	for(int d = 0; d < dim; d++){		
		std::vector<double> checkedEdges = {0};
		//For each edge
		for(auto edge : local_edges){
			
			//Get the weights (increasing order)
			double epsilon = edge.first;
			//Check if we've already 
			if(std::find(checkedEdges.begin(), checkedEdges.end(), epsilon) == checkedEdges.end()){
				
				std::cout << "EPSILON: " << epsilon << std::endl;
				//Get the reduced boundary matrix
				std::pair<std::vector<std::vector<unsigned>>, int> boundary = boundaryMatrix(nSimplices(epsilon,d+1,local_edges),nSimplices(epsilon,d,local_edges));
				if(boundary.second > 0){
					
					std::cout << "\tReduced boundary matrix @ " << epsilon << " -> " << boundary.first.size() ;
					std::cout << " x " <<boundary.first[1].size() << std::endl;
					
					std::cout << std::endl << "REDUCED BOUNDARY" << std::endl;
					for(unsigned z = 0; z < boundary.first.size(); z++)
						ut.print1DVector(boundary.first[z]);
					std::cout << std::endl << std::endl;
					extractBoundaries(local_edges,boundary.first, boundary.second);
				
					std::cout << "\tDim: " << d << "\tRANKS: ";
					for(auto z : ranks)
						std::cout << z << " ";
					std::cout << std::endl;
				}
		
				checkedEdges.push_back(epsilon);
				allBoundaries.push_back(boundary.first);
			}
		}
		checkedEdges={0};
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

