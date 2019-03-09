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
//		generate boundary matrix, determine constituent boundaries
//		edges that make up the boundary
//
//		For each dimension, 
//			Generate boundary matrix at each weight
//				Then, reduce the boundary matrix to RREF
//				Extract boundaries from RREF
//			Generate barcodes (lifespans) of each component
//				WITH boundaries that form the barcode
//
pipePacket boundaryPipe::runPipe(pipePacket inData){
	std::vector<std::vector<unsigned>> allBoundaries;
	
	
	struct bettiDef_t{
		double epsilon;
		int dim;
		int betti;
	};
	
	std::vector<bettiDef_t> bettis;
	
	
	std::vector<int> bettiNumbers;
	std::vector<float> lifeSpans[dim];
	
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges;
	
	//Retrieve
	auto local_weights = inData.weights;
	std::string barcodes;
	for(int d2 = 0; d2 <= dim; d2++){
		bettiNumbers.push_back(0);
	}
	
	std::string bettiOutput[] = {"Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n","Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n","Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n"};
	
	double epsilon = 0;
	
	//For each edge
	for(auto eps : local_weights){
		//Reload the buffers with the current edges
		edges = inData.workData.complex->getAllEdges(eps);
		
		//Get the weights (increasing order)
		//Check if we've already processed or not
		if(epsilon != eps){
			epsilon = eps;
			auto last_rank_nul = std::make_pair(0,0);
			
			//Iterate through each dimension to build boundary matrix
			for(int d = dim; d >= 0; d--){
				
				//Get the reduced boundary matrix
				std::pair<std::vector<std::vector<unsigned>>,std::pair<int,int>>  bound_rank_nul;
				
				//if(d == 0)
				//	bound_rank_nul = boundaryMatrix({}, edges[d]);
				//else
				//	bound_rank_nul = boundaryMatrix(edges[d-1], edges[d]);
				
				if(d > 1){
					std::cout <<"Extracting boundaries..." << std::endl;
					/*std::vector<std::vector<unsigned>> z = extractBoundaries(edges[d-1], bound_rank_nul.first, bound_rank_nul.second.second);
					
					std::cout << "\n\n______________BOUNDARIES (" << std::to_string(z.size()) << ")_______________" << std::endl;
					
					for(auto a : z){
						ut.print1DVector(a);
					}
					
					
					for(auto bound : z){
						std::cout << "bound\t";
						for(auto curBound : allBoundaries){
							std::cout << "curBound\t";
							if(ut.setIntersect(bound, curBound, true).size() > bound.size())
								allBoundaries.push_back(bound);
						}
					}*/
				}
				
				auto rank_nul = bound_rank_nul.second;
											
				if(bettiNumbers[d] != (rank_nul.second- last_rank_nul.first)){
					bettiNumbers[d] = (rank_nul.second- last_rank_nul.first);
				}
				
				if(d != dim)
					bettis.push_back(bettiDef_t {epsilon, d, rank_nul.second - last_rank_nul.first});
				if(debug){
					bettiOutput[d] += std::to_string(epsilon) + "\t" + std::to_string(d) + "\t" + std::to_string(d-1) + "\t" + std::to_string(rank_nul.second)/* - last_rank_nul.first)*/ + "\t" + std::to_string(rank_nul.first) + "\t" + std::to_string(rank_nul.second) + "\t" + std::to_string(last_rank_nul.first) + "\t" + std::to_string(last_rank_nul.second) + "\n";
				}
				last_rank_nul = rank_nul;
				
			}
		}

	}
	
	std::string output = "Dim,Birth,Death\n";
	
	for(int i = 0; i < dim; i ++){
		bettiDef_t lastBetti = {0.0,i,0};
		for(int j = bettis.size(); j >= 0; j--){
			int lastBirth;
			
			if(bettis[j].dim == i){
				
				if(bettis[j].betti < lastBetti.betti){
					for(int k = 0; k < (lastBetti.betti - bettis[j].betti); k++)
						output += std::to_string(i) + "," + std::to_string(lastBetti.epsilon) + "," + std::to_string(bettis[j].epsilon) + "\n";
					lastBetti.betti = bettis[j].betti;
				}
				if(bettis[j].betti > lastBetti.betti){
					lastBetti.epsilon = bettis[j].epsilon;
					lastBetti.betti = bettis[j].betti;
				}
			}
			
			if(j == 0){
				while(lastBetti.betti > 0){
					output += std::to_string(i) + "," + std::to_string(lastBetti.epsilon) + "," + std::to_string(maxEpsilon) + "\n";
					lastBetti.betti--;
				}
			}
		}
	}
	
	if(debug){
		std::cout << "\n\n______________RESULTS_______________" << std::endl;
		std::cout << bettiOutput[0] << std::endl << std::endl;
		std::cout << bettiOutput[1] << std::endl << std::endl;
		std::cout << bettiOutput[2] << std::endl;
		std::cout << std::endl << output << std::endl;
	}
	
	std::cout << "\n\n______________BOUNDARIES_______________" << std::endl;
	for(auto a : allBoundaries){
		ut.print1DVector(a);
	}
	
	inData.bettiOutput = output;
	
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
	
	return true;
}

