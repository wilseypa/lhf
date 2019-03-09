/*
 * bettiPipe hpp + cpp extend the basePipe class for calculating the 
 * betti numbers from a distance matrix
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

// Check if a face is a subset of a simplex
int bettiPipe::checkFace(std::set<unsigned> face, std::set<unsigned> simplex){
	
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

// Reduce the binary boundary matrix
std::pair<std::queue<unsigned>,std::pair<int,int>> bettiPipe::reduceBoundaryMatrixRev(std::vector<std::vector<unsigned>> boundaryMatrix){
	std::vector<std::vector<unsigned>> ret;
	std::queue<unsigned> pivots;
	int rank = 0;
	int invNul = boundaryMatrix.size();
	if(boundaryMatrix.size() == 0 || boundaryMatrix[0].size() <= 0)
		return std::make_pair(pivots, std::make_pair(0,0));
	
	//Step through each column and search for a 1 in that column
	
	for(unsigned p = 0; p < boundaryMatrix.size(); p++){
		//Step through each vector
		for(unsigned n = 0; n < boundaryMatrix[0].size(); n++){
			//If the vector has a 1 in the target column
			if(boundaryMatrix[p][n] == 1){
				rank += 1;
				pivots.push(p);
				
				//Add the vector to our returned matrix
				int tempRow = n;
					
				//Remove our vector from the boundaryMatrix and continue processing next column
				
				//Iterate through remaining vectors and XOR with our vector or the clearRow
				//j += 1;
				while(n < boundaryMatrix[0].size()){
					if(boundaryMatrix[p][n] == 1){
						//XOR the remaining rows
						for(unsigned d = 0; d < boundaryMatrix.size(); d++)
							boundaryMatrix[d][n] = boundaryMatrix[d][n] ^ boundaryMatrix[d][tempRow];
					}
					n += 1;
				}			
				
				for(auto bm : boundaryMatrix)
					bm.erase(bm.begin() + tempRow);
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
	
	return std::make_pair(pivots, std::make_pair(rank, (invNul - rank)));
}


// Get the boundary matrix from the edges
std::pair<std::queue<unsigned>, std::pair<int,int>> bettiPipe::getRank(std::vector<std::vector<unsigned>> nChain, std::vector<std::vector<unsigned>> pChain, std::queue<unsigned> pivots){
	std::queue<unsigned> temp;
	if (nChain.size() == 0)
		return std::make_pair(temp, std::make_pair(0,pChain.size()));
	
	if(twist == "true"){
		std::vector<std::vector<unsigned>> boundary;
		//Clear out the pChain columns faces of zChain from the boundary matrix 
		//	See [Chen, Kerber 2011] for details
		//	
		//	When the boundary matrix is created, clear out columns
		//	that are a coface of the higher level simplice list
		//		This requires the list of simplices be passed down
		//		BUT it reduces the complexity of reduction
		//			(may require restructuring of boundary matrix (pivoted))
		
		
		unsigned curPivot = -1;
		if(!pivots.empty()){
			curPivot = pivots.front();
			pivots.pop();
		}
			
		int clearedCount = 0;
		//Create the columns (pChain)
		for(unsigned i = 0; i < pChain.size(); i++){
			std::vector<unsigned> tempVector(nChain.size(), 0);
			//Check if this is a
				
			if(i == curPivot){
					clearedCount++;
					if(!pivots.empty()){
						curPivot = pivots.front();
						pivots.pop();
					}
			}			
			else {
				//Create the columns (pChain)
				for(unsigned j = 0; j < nChain.size(); j++){
					tempVector[j] = (checkFace(nChain[j], pChain[i]));
				}
			}
			
			boundary.push_back(tempVector);
		}	
		if(debug)
			std::cout << "\tCleared Count: " << clearedCount << std::endl;
		
		//Reduce the boundary matrix
		return reduceBoundaryMatrixRev(boundary);
	}

	std::vector<std::vector<unsigned>> boundary;

	
	
	//Reduce the boundary matrix
	return std::make_pair(temp, reduceBoundaryMatrix(boundary));
}

std::pair<std::queue<unsigned>, std::pair<int,int>> bettiPipe::removeSimplices(std::vector<std::pair<std::set<unsigned>,double>> &dimEdges, std::vector<std::vector<unsigned>> &dimBoundary, double epsilon, std::queue<unsigned> &pivots, std::pair<int, int> &lastRankNull){
	//Determine if there are simplices to remove from this dimension at the new epsilon
	//		If there are, update the edges and boundary matrix for the dimension; return true (as in recompute)
	//		Otherwise return false; no update is needed
	std::queue<unsigned> retPivots;
	std::pair<std::queue<unsigned>, std::pair<int,int>> ret;
	
	//This is the first iteration, don't want to remove anything and not worth iterating
	if(epsilon == maxEpsilon){
		if(twist == "true"){
			return reduceBoundaryMatrixRev(dimBoundary);
		} else {
			//Not using pivots
			return std::make_pair(retPivots, reduceBoundaryMatrix(dimBoundary));
		}
	}
	
	//Store the edges to remove to prevent inconsistencies when indexing
	std::vector<std::pair<std::set<unsigned>,double>> remEdges;
	
	//Track the number of removed items for consistency when indexing
	int removedItems = 0;
	
	//Iterate through each of the edges, looking for anything cut off by epsilon
	//		these entries should be sorted desc by weight
	for(int i = 0; i < dimEdges.size(); i++){
		
		//If we've reached edges that remain
		if(dimEdges[i].second <= epsilon)
			break;
			
		remEdges.push_back(dimEdges[i]);
		removedItems++;		
	}
	
	//Remove the edges from the boundary matrix
	if(twist == "true"){
		if(removedItems > 0){
			//When twisted, remove the rows (from start to # of removed items)
			dimBoundary.erase(dimBoundary.begin(), dimBoundary.begin() + removedItems);
			
			ret = reduceBoundaryMatrixRev(dimBoundary);
		} else {
			ret = std::make_pair(retPivots, lastRankNull);
		}
			
	} else {
		if(removedItems > 0){
			//Otherwise remove the columns (from start to # of removed items)
			for (int i = 0; i < dimBoundary.size(); i++)
				dimBoundary[i].erase(dimBoundary[i].begin(), dimBoundary[i].begin() + removedItems);
				
			ret = std::make_pair(retPivots, reduceBoundaryMatrix(dimBoundary));
		} else {
			
			ret = std::make_pair(retPivots, lastRankNull);
		}
			
	}
	//Remove the edges from the edge list
	dimEdges.erase(dimEdges.begin(), dimEdges.begin() + removedItems);
	
	return ret;
	
}

std::vector<std::vector<std::vector<unsigned>>> bettiPipe::createBoundaryMatrix(std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges){
	std::vector<std::vector<std::vector<unsigned>>> boundary;
	
	//Iterate through each dimensional graph...
	for(int d = edges.size()-1; d >= 0; d--){
		std::vector<std::vector<unsigned>> tempBoundary;
		
		//Setup p-chains and n-chains
		std::vector<std::pair<std::set<unsigned>,double>> nChain;
		std::vector<std::pair<std::set<unsigned>,double>> pChain = edges[d];
		
		
		if(d == 0)
			nChain = {};
		else
			nChain = edges[d-1];
			
			
		//Decide if we're going to twist or not; changes orientation of the boundary matrix
		if(twist == "true"){
		
			//Create the columns (pChain) PIVOTED!
			for(unsigned i = 0; i < pChain.size(); i++){
				std::vector<unsigned> tempVector(nChain.size(), 0);
				
				//Create the rows (nChain)
				for(unsigned j = 0; j < nChain.size(); j++){
					tempVector[j] = (checkFace(nChain[j].first, pChain[i].first));
				}
				tempBoundary.push_back(tempVector);	
			}
		} else {
			//Create the rows (nChain)
			for(unsigned i = 0; i < nChain.size(); i++){
			std::vector<unsigned> bDim;
			
				//Create the columns (pChain)
				for(unsigned j = 0; j < pChain.size(); j++){
					bDim.push_back(checkFace(nChain[i].first, pChain[j].first));
				}
				tempBoundary.push_back(bDim);
			}	
		}

		
		boundary.insert(boundary.begin(),tempBoundary);
		
	}

	return boundary;
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
	
	
	std::vector<int> bettiNumbers (dim,0);
	std::vector<float> lifeSpans[dim];
	
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges = inData.workData.complex->getAllEdges(maxEpsilon);
	std::vector<std::vector<std::vector<unsigned>>> dimBoundaries = createBoundaryMatrix(edges);
	std::queue<unsigned> dimPivots[dim+1];
	
	std::pair<int, int> dimRankNull[dim+2];
	std::fill_n(dimRankNull, dim+1, std::make_pair(0,0));
	
	//Retrieve
	auto local_weights = inData.weights;
	std::string barcodes;
	
	std::string bettiOutput[dim+1];
	std::fill_n(bettiOutput, dim+1, "Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n");
	
	double epsilon = 0;
	
	//For each edge
	for(auto eps : local_weights){
		//Reload the buffers with the current edges
		//	note - this is no longer controlled from the complex structure;
		//
		//	Instead, we should look at the previous boundaries (dn, ...., d0) 
		//		for changes; if the weight is < epsilon for a dimension we need
		//		to recompute the boundaries (but only remove a feature; don't 
		//		recreate the boundary matrix from scratch)
		//
		//	To do this, store each boundary matrix locally; pass by reference
		//		to the rank function if the boundary matrix has changed
		
		//Get the weights (increasing order)
		//Check if we've already processed or not
		if(epsilon != eps){
			epsilon = eps;
			
			auto last_rank_nul = std::make_pair(0,0);
			std::queue<unsigned> pivots;
			//Iterate through each dimension to build boundary matrix
			for(int d = dim; d >= 0; d--){
				//If this is true, the rank/nullity of the dimension needs
				//	to be recalculated; the boundary matrix has changed!
				std::pair<std::queue<unsigned>, std::pair<int,int>> retVal;
				
				if(d == 0)
					retVal.second = std::make_pair(0, edges[d].size());
				else
					retVal = removeSimplices(edges[d], dimBoundaries[d], epsilon, dimPivots[d+1], dimRankNull[d]);
				
				
				dimRankNull[d] = retVal.second;
				dimPivots[d] = retVal.first;
				
				if(d != dim)
					bettis.push_back(bettiDef_t {epsilon, d, dimRankNull[d].second - dimRankNull[d+1].first});
				if(debug){
					bettiOutput[d] += std::to_string(epsilon) + "\t" + std::to_string(d) + "\t" + std::to_string(d-1) + "\t" + std::to_string(dimRankNull[d].second - dimRankNull[d+1].first)/* - last_rank_nul.first)*/ + "\t" + std::to_string(dimRankNull[d].first) + "\t" + std::to_string(dimRankNull[d].second) + "\t" + std::to_string(dimRankNull[d+1].first) + "\t" + std::to_string(dimRankNull[d+1].second) + "\n";
				}	
				
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
		for(auto z : bettiOutput){
			std::cout << z << std::endl << std::endl;
		}
		std::cout << std::endl << output << std::endl;
	}
	
	inData.bettiOutput = output;
	
	return inData;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void bettiPipe::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_bettis_output.csv");
	
	file << inData.bettiOutput;
	
	file.close();
	return;
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
	
	pipe = configMap.find("twist");
	if(pipe != configMap.end())
		twist = configMap["twist"];
	else return false;
	
	return true;
}

