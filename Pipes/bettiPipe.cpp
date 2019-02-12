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


// basePipe constructor
bettiPipe::bettiPipe(){
	pipeType = "Betti";
	return;
}

//Filter and return simplices of a specified dimension
std::vector<std::vector<unsigned>> bettiPipe::nSimplices(double epsilon, unsigned n, std::vector<std::pair<double,std::vector<unsigned>>> complex){
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
int bettiPipe::checkFace(std::vector<unsigned> face, std::vector<unsigned> simplex){
	
	if(simplex.size() == 0)
		return 1;
	else if(ut.symmetricDiff(face,simplex,false).size() == 1){
		return 1;
	}
	else
		return 0;
}

// Reduce the binary boundary matrix
std::pair<int,int> bettiPipe::reduceBoundaryMatrix(std::vector<std::vector<unsigned>> boundaryMatrix){
	std::vector<std::vector<unsigned>> ret;
	int rank = 0;
	
	if(boundaryMatrix.size() <= 0)
		return std::make_pair(0,0);
		
	//Step through each column and search for a 1 in that column
	for(unsigned i = 0; i < boundaryMatrix[0].size(); i++){
		//Step through each vector
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If the vector has a 1 in the target column
			if(boundaryMatrix[j][i] == 1){
				rank += 1;
				
				//Add the vector to our returned matrix
				auto tempRow = boundaryMatrix[j];
					
				//Remove our vector from the boundaryMatrix and continue processing next column
				boundaryMatrix.erase(boundaryMatrix.begin() + j);
				
				//Iterate through remaining vectors and XOR with our vector or the clearRow
				//j += 1;
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j][i] == 1){
						//XOR the remaining rows
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] ^ tempRow[d];
					}
					j += 1;
				}			
				
				/*for(int a = 0; a < ret.size(); a++){
					//xor the returned rows
					if (ret[a][i] == 1){
						for(unsigned d = 0; d < ret[a].size(); d++)
							ret[a][d] = (ret[a][d] ^ tempRow[d]);
					}					
				}
				ret.push_back(tempRow);*/
			}
		}
	}

	if(boundaryMatrix.size() > 0){
		for (auto a : boundaryMatrix){
			ret.push_back(a);
		}
	}
	
	return std::make_pair(rank, (ret[0].size() - rank));
}

// Get the boundary matrix from the edges
std::pair<int,int> bettiPipe::getRank(std::vector<std::vector<unsigned>> nChain, std::vector<std::vector<unsigned>> pChain){
	std::vector<std::vector<unsigned>> boundary;	

	if (nChain.size() == 0)
		return std::make_pair(0,pChain.size());
	
	//Create the rows (nChain)
	for(unsigned i = 0; i < nChain.size(); i++){
		std::vector<unsigned> bDim;
		
		//Create the columns (pChain)
		for(unsigned j = 0; j < pChain.size(); j++){
			bDim.push_back(checkFace(nChain[i], pChain[j]));
		}
		boundary.push_back(bDim);
	}	
	
	//Reduce the boundary matrix
	return reduceBoundaryMatrix(boundary);
}


// runPipe -> Run the configured functions of this pipeline segment
//
//	bettiPipe: For computing the betti numbers from simplicial complex:
//		1. Compute betti numbers from ranks at each weight of complex
//			a. Point to sorted simplical complex and get next weight (epsilon)
//			b. Retrieve all points less than or equal to weight
//			c. Compute betti numbers up to max dimension
//			d. Check if betti numbers have changed since last iteration;
//				Record changes as birth times/death times
//		
//
pipePacket bettiPipe::runPipe(pipePacket inData){
	
	struct bettiDef_t{
		double epsilon;
		int dim;
		int betti;
	};
	
	std::vector<bettiDef_t> bettis;
	
	std::vector<int> bettiNumbers;
	std::vector<float> lifeSpans[dim];
	
	std::vector<std::vector<std::vector<unsigned>>> edges;
	
	//Retrieve
	auto local_weights = inData.workData.complex->weights;
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
				std::pair<int, int> rank_nul;
				
				if(d == 0)
					rank_nul = getRank({}, edges[d]);
				else
					rank_nul = getRank(edges[d-1], edges[d]);
							
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
	
	if(debug){
		std::cout << "\n\n______________RESULTS_______________" << std::endl;
		std::cout << bettiOutput[0] << std::endl << std::endl;
		std::cout << bettiOutput[1] << std::endl << std::endl;
		std::cout << bettiOutput[2] << std::endl;
	}
	
	std::cout << "\n\n______________RESULTS_______________" << std::endl;
	
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
	
	std::cout << output << std::endl;
	
	
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
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	return true;
}

