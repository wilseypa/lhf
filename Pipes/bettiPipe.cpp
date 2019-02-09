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
		
		//for(auto i : v.second)
		//	std::cout << i << " ";
		//std::cout << "--> " << v.second.size() << " " << v.first << " " << epsilon << std::endl << std::endl;
		
		
		if(v.second.size() == n){
			if(v.first <= epsilon)
				ret.push_back(v.second);
		}
	}
	
	/*std::cout << "nSimplices : " << n << std::endl;
	for(auto i : ret){
		for(auto p : i)
			std::cout << p << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;*/
	
	return ret;
}

// Check if a face is a subset of a simplex
int bettiPipe::checkFace(std::vector<unsigned> face, std::vector<unsigned> simplex){
	/*std::cout << std::endl << "________CheckFace____________" << std::endl;
	ut.print1DVector(face);
	ut.print1DVector(simplex);*/
	//std::cout << "SIZE: " << simplex.size() << std::endl;
	//ut.print1DVector(ut.intersect(face,simplex,false).second);

	
	if(simplex.size() == 0)
		return 1;
	else if(ut.symmetricDiff(face,simplex,false).size() == 1){
		//ut.print1DVector(face);
		//ut.print1DVector(simplex);
		return 1;
	}
	else
		return 0;
}

// Reduce the binary boundary matrix
std::pair<int,int> bettiPipe::reduceBoundaryMatrix(std::vector<std::vector<unsigned>> boundaryMatrix, int nChainSize){
	std::vector<std::vector<unsigned>> ret;
	int rank = 0;
	if(boundaryMatrix.size() <= 0)
		return std::make_pair(0,0);
		
	if(nChainSize > 100000){
		nChainSize = 1;
	}
	
	/*std::cout << std::endl;
	if(boundaryMatrix.size() > 0){
		for (auto a : boundaryMatrix){
			ut.print1DVector(a);
		}
	} 
	std::cout << "-------------------------------------------" << std::endl;
		*/
	//Step through each column and search for a 1 in that column
	for(unsigned i = 0; i < boundaryMatrix[0].size(); i++){
		bool found = false;
		//Step through each vector
		for(unsigned j = 0; j < boundaryMatrix.size(); j++){
			
			//If the vector has a 1 in the target column
			if(boundaryMatrix[j][i] == 1){
				rank += 1;
				found = true;
				
				//Add the vector to our returned matrix
				//Along with a vector to clear out the row
				auto tempRow = boundaryMatrix[j];
				//auto clearRow = boundaryMatrix[j];
				//clearRow[i] = 0;
				
				//Set up the temp row for both operations:
				
				//for(unsigned d = i+1; d < boundaryMatrix[0].size(); d++)
				//	clearRow[d] = tempRow[d];
					
				//Remove our vector from the boundaryMatrix and continue processing next column
				boundaryMatrix.erase(boundaryMatrix.begin() + j);
				
				//Iterate through remaining vectors and XOR with our vector or the clearRow
				//j += 1;
				while(j < boundaryMatrix.size()){
					if(boundaryMatrix[j][i] == 1){
						
						//XOR the two rows
						for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							boundaryMatrix[j][d] = boundaryMatrix[j][d] != tempRow[d];
					}
					//else{
						//for(unsigned d = 0; d < boundaryMatrix[0].size(); d++)
							//boundaryMatrix[j][d] = boundaryMatrix[j][d] != clearRow[d];
					//}
							
					j += 1;
				}			
				
				for(int a = 0; a < ret.size(); a++){
					//xor the returned rows
					if (ret[a][i] == 1){
						std::cout << "_________" << std::endl;
						ut.print1DVector(tempRow);
						ut.print1DVector(ret[a]);
						for(unsigned d = 0; d < ret[a].size(); d++)
							ret[a][d] = (ret[a][d] != tempRow[d]);
						ut.print1DVector(ret[a]);
					}					
				}
				
				std::cout << "_________" << std::endl;
				ret.push_back(tempRow);
			}
		}
	}

	if(boundaryMatrix.size() > 0){
		for (auto a : boundaryMatrix){
			ret.push_back(a);
		}
	}
	
	for(auto z: ret){
		ut.print1DVector(z);
	}
	
	std::cout << "ChainSize: " << nChainSize << std::endl;
	std::cout << "RetSize: " << ret.size() << std::endl;
	std::cout << "Rank: " << rank << std::endl;
	if(ret[0].size() > 1){
		
		std::cout << "Nullity: " << (ret[0].size() - rank) << std::endl;
	}
	ranks.push_back(rank);
	nRanks.push_back(ret.size() - rank);
	
	if(ret.size() > 1){
		return std::make_pair(rank, (ret.size() - rank));
	}
	else
		return std::make_pair(rank, ret[0].size() - rank);
}

// Get the boundary matrix from the edges
std::pair<int,int> bettiPipe::getRank(std::vector<std::vector<unsigned>> nChain, std::vector<std::vector<unsigned>> pChain, int nCount){
	std::vector<std::vector<unsigned>> ret;
	std::vector<std::vector<unsigned>> boundary;	
	
	
	//Create the boundary matrix from chains
	std::vector<unsigned> bDim;
	
	std::cout << "pre: " << nChain.size() << " " << pChain.size() << std::endl;
	
	if (nChain.size() == 0)
		return std::make_pair(0,pChain.size());
	
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
		/*
		std::cout << "_______pChain______" << std::endl;
		for(auto i : pChain){
			for(auto p : i)
				std::cout << p << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;

		*/
		//Create the rows (pChain)
		for(unsigned i = 0; i < pChain.size(); i++){
			//Create the columns (nChain)
			for(unsigned j = 0; j < nChain.size(); j++){
				bDim.push_back(checkFace(pChain[i], nChain[j]));
			}
			boundary.push_back(bDim);
			bDim.clear();
		}	
		std::cout << std::endl;
	}
		
	std::cout << "\tFinished boundary matrix -> " << boundary.size() << " x " <<boundary[1].size() << std::endl << std::endl;
	
	
	std::cout << std::endl << "ORIGINAL BOUNDARY" << std::endl;
	for(unsigned d = 0; d < boundary.size(); d++)
		ut.print1DVector(boundary[d]);
	std::cout << std::endl << std::endl;
	
	//Reduce the boundary matrix
	return reduceBoundaryMatrix(boundary, nChain.size());
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
	std::vector<std::vector<std::vector<unsigned>>> allBoundaries;
	
	std::vector<int> bettiNumbers;
	std::string bettiOutput[] = {"Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n","Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n","Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n"};
	
	//Retrieve
	auto local_edges = inData.workData.complex->getEdges(0,0);
	std::string barcodes;
	for(int d2 = 0; d2 <= dim; d2++){
		bettiNumbers.push_back(0);
	}
	
	double epsilon = 0;
	
	//For each edge
	for(auto edge : local_edges){
		
		//Get the weights (increasing order)
		//Check if we've already processed or not
		if(epsilon != edge.first){
			epsilon = edge.first;
			
			auto last_rank_nul = std::make_pair(0,0);
			
			//Iterate through each dimension to build boundary matrix
			for(int d = dim; d >= 0; d--){
				
				std::cout << "Test-" << d << "================================>>" << std::endl;
				//Get the reduced boundary matrix
				auto rank_nul = getRank(nSimplices(epsilon,d+1,local_edges),nSimplices(epsilon,d,local_edges), nSimplices(epsilon,d+2,local_edges).size());
				
				//If the rank is greater than 0 (i.e. a feature exists in the boundary matrix)
				std::cout << "\t" << rank_nul.first << "," << rank_nul.second << "\td" << d << "\t" ;
				
				std::cout << "\tBetti " << d << ": " << (rank_nul.second - last_rank_nul.first) << std::endl;
				std::cout << "\tEpsilon " << epsilon << std::endl;
				std::cout << "LastNull: " << last_rank_nul.first << std::endl;
				std::cout << "curNull: " << rank_nul.second << std::endl;
				
				
				//if(bettiNumbers[d] != (rank_nul.second- last_rank_nul.first)){
					std::cout << "Test" << d << std::endl;
					bettiOutput[d] += std::to_string(epsilon) + "\t" + std::to_string(d) + "\t" + std::to_string(d-1) + "\t" + std::to_string(rank_nul.second)/* - last_rank_nul.first)*/ + "\t" + std::to_string(rank_nul.first) + "\t" + std::to_string(rank_nul.second) + "\t" + std::to_string(last_rank_nul.first) + "\t" + std::to_string(last_rank_nul.second) + "\n";
					
					bettiNumbers[d] = (rank_nul.second- last_rank_nul.first);
				//}
				
				last_rank_nul = rank_nul;
				
			}
			std::cout << "\t| " << epsilon << std::endl;
		}

	}
	
	std::cout << "\n\n______________RESULTS_______________" << std::endl;
	std::cout << bettiOutput[0] << std::endl << std::endl;
	std::cout << bettiOutput[1] << std::endl << std::endl;
	std::cout << bettiOutput[2] << std::endl;

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

