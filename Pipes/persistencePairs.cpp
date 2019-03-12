/*
 * persistencePairs hpp + cpp extend the basePipe class for calculating the 
 * persistence pairs numbers from a complex
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
#include "persistencePairs.hpp"


// basePipe constructor
persistencePairs::persistencePairs(){
	pipeType = "PersistencePairs";
	return;
}

bool sortBySecondAsc(const std::pair<std::set<unsigned>, double> &a, const std::pair<std::set<unsigned>, double> &b){
	return (a.second < b.second);
}
	
//Filter and return simplices of a specified dimension
std::vector<std::vector<unsigned>> persistencePairs::nSimplices(double epsilon, unsigned n, std::vector<std::pair<double,std::vector<unsigned>>> complex){
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
int persistencePairs::checkFace(std::vector<unsigned> face, std::vector<unsigned> simplex){
	
	if(simplex.size() == 0)
		return 1;
	else if(ut.symmetricDiff(face,simplex,false).size() == 1){
		return 1;
	}
	else
		return 0;
}

// Check if a face is a subset of a simplex
int persistencePairs::checkFace(std::set<unsigned> face, std::set<unsigned> simplex){
	
	if(simplex.size() == 0)
		return 1;
	else if(ut.symmetricDiff(face,simplex,false).size() == 1){
		return 1;
	}
	else
		return 0;
}

std::vector<unsigned> persistencePairs::getRankNull(std::vector<std::vector<unsigned>> boundaryMatrix, unsigned offset){
	
	//Perform column echelon reduction; basically the inverse of the RREF
	//	return offset index of pivots (ranks)
	
	std::vector<unsigned> ret;
	
	//Step through each row and search for a 1 in that column
	for(unsigned i = 0; i < boundaryMatrix.size(); i++){
		
		//Step through each column
		for(unsigned j = 0; j < boundaryMatrix[0].size(); j++){
			
			//If the vector has a 1 in the target column
			if(boundaryMatrix[i][j] == 1){
				ret.push_back(j + offset);
	
				//Create a temp matrix for XOR
				std::vector<unsigned> tempColumn;
				
				//Zero out our column and continue processing next column to clear any ones
				for(int z = 0; z < boundaryMatrix.size(); z++){
					tempColumn.push_back(boundaryMatrix[z][j]);
					boundaryMatrix[z][j] = 0;
				}
				
				//Iterate through remaining vectors and XOR with our vector or the clearRow
				//j += 1;
				while(j < boundaryMatrix[0].size()){
					if(boundaryMatrix[i][j] == 1){
						//XOR the remaining rows
						for(unsigned d = 0; d < boundaryMatrix.size(); d++)
							boundaryMatrix[d][j] = boundaryMatrix[d][j] ^ tempColumn[d];
					}
					j += 1;
				}
			}
		}
	}
	
	return ret;
}

std::vector<unsigned> persistencePairs::createBoundaryMatrix(std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges){
	std::vector<unsigned> pivots;
	
	int totalEdges = 0;
	for(auto a : edges)
		totalEdges += a.size();
	
	int lastCount = 0;
	
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
			
		//Create the rows (nChain)
		for(unsigned i = 0; i < nChain.size(); i++){
		std::vector<unsigned> bDim;
		
			//Create the columns (pChain)
			for(unsigned j = 0; j < pChain.size(); j++){
				bDim.push_back(checkFace(nChain[i].first, pChain[j].first));
			}
			tempBoundary.push_back(bDim);
		}	
		
		//Get the offset pivots of the boundary matrix in column echelon form
		auto tempPivots = getRankNull(tempBoundary, totalEdges - tempBoundary.size() - lastCount - 1);	
		
		lastCount = totalEdges - tempBoundary.size();
		
		for(auto a : tempPivots)
			pivots.push_back(a);
		
		
	}

	
	ut.print1DVector(pivots);

	return pivots;
}


std::vector<std::vector<std::pair<double,double>>> persistencePairs::computeIntervals(std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges){
	std::vector<std::vector<std::pair<double,double>>> ret;
	std::vector<std::pair<double,double>> temp;
	for(int i = 0; i < dim; i++){
		ret.push_back(temp);
	}
		
	//def ComputeIntervals(K):
	//for k=0 to dim(K) Lk = 0; {
	//	for j=0 to m-1 {
	//		d = RemovePivotrows(sj)
	//		if(d==0) Mark sj
	//		else {
	//			i = maxindexd; k = dim si;
	//			store j and d in T[i];
	//		}
	//	}
	//	for j=0 to m-1 {
	//		if sj is marked and T[j] is empty{
	//			k = dim sj; Lk = Lk U {(deg sj, inf)}
	//		}
	//	}
	//}
	
	//In basic terms
	//		Iterate through total-ordered simplices, (by dim, weight)
	//			Call RemovePivotRows and retrieve simplex or null
	//				if null, mark the simplex (in T[i]) and continue
	//				if there is a ret, find that simplex in the T-array (T[i])
	//					Generate a persistence interval (weight of T[i], weight of T[j])
	//					Store j and d in T[i];
	//		Afterwards, find all marked rows, j, with no entry for T[j]
	//			Generate persistence interval for these, (weight of T[j], inf)
	//		Success!
	
	//Generate our boundary matrix and retrieve pivots
	std::cout << "Creating boundary matrix... " << std::endl;
	std::vector<unsigned> kPivots = createBoundaryMatrix(edges);
	
	std::string bettis = "curIndex,index,dim,birth,death\n";
	
	//Flatten the edges into a single array
	std::vector<std::set<unsigned>> kSimplices;
	std::vector<double> kWeights;
	for(auto a : edges){		
		for(auto z : a){
			kSimplices.push_back(z.first);
			kWeights.push_back(z.second);
		}
	}
	
	//Create our T-array
	std::vector<tArrayEntry_t> tArray;
	
	std::cout << "Compute intervals..." << std::endl;
	
	for(int curIndex = 0; curIndex < kSimplices.size(); curIndex++){
		int curDim = kSimplices[curIndex].size();
		tArray.push_back(tArrayEntry_t());
		
		if(curDim == 1){
			tArray[curIndex].birth = kWeights[curIndex];
			tArray[curIndex].marked = true;
		}
		else if(std::find(kPivots.begin(), kPivots.end(), curIndex) != kPivots.end()){
			std::cout << "Found " << curIndex << " in d = " << curDim - 1 << std::endl;
			tArray[curIndex].birth = kWeights[curIndex];
			tArray[curIndex].marked = true;
		} else {
			std::set<unsigned> simp = kSimplices[curIndex];
			
			int iter = tArray.size();
			for(; iter >= 0; iter--){
				if(tArray[iter].marked && tArray[iter].death < 0 && simp.size() != kSimplices[iter].size() && ut.setIntersect(kSimplices[iter], simp, false).size() == simp.size() - 1)
					break;
			}
			
			std::cout << "Subtree..." << std::endl;
			
			//Get simp_i
			//std::vector<std::set<unsigned>>::iterator kIndex = std::find(kSimplices.begin(), kSimplices.end(), simp);
			
			if(tArray[iter].death < 0 ){
				tArray[iter].death = kWeights[curIndex];
				
				bettis += std::to_string(curIndex) + "," + std::to_string(iter) + "," + std::to_string(simp.size() - 2) + "," + std::to_string(tArray[iter].birth) + "," + std::to_string(kWeights[curIndex]) + "\n";
			}
			
		}
	}
	
	std::cout << "Br" << std::endl;
	std::cout << bettis << std::endl;
	
	std::cout << "T\tWeights\tti\tbirth\tdeath\tmarked\tsize" << std::endl;
	for(int t = 0; t < kSimplices.size(); t++){
		std::cout << t << "\t" << kWeights[t] << "\t" << tArray[t].ti << "\t" << tArray[t].birth << "\t" << tArray[t].death << "\t" << tArray[t].marked << "\t" << tArray.size() << "\t";
		ut.print1DVector(kSimplices[t]);
		if(tArray[t].marked && tArray[t].ti == -1){
			ret[tArray[t].simplex.size()].push_back(std::make_pair(kWeights[t], maxEpsilon));
		}
	}
			
		
	
	//for k=0 to dim(K) Lk = 0;
	/*for(auto simp : kSimplices){													//	for j=0 to m-1 {
		std::set<unsigned> d = removePivotRows(simp, dimBoundaries, tArray);								//		d = RemovePivotrows(sj)
		if(d.size() == 0)
			tArray[curIndex].marked = true;											//		if(d==0) Mark sj
		else {
			unsigned i = *d.rbegin();												//			i = maxindexd; k = dim si;
			auto k = kSimplices[i].first.size();									// !! need to get kSimplices[i].size in a diff manner
			tArray[i].maxIndex = curIndex;
			tArray[i].simplex = d;													//			store j and d in T[i];
			
			intervals[simp.first.size()].push_back(std::make_pair(kSimplices[i].second, kSimplices[curIndex].second));																//NOTE: Degree is the weight of the simplex
		}
		curIndex++;
	}
	
	curIndex = 0;
	for(auto simp : kSimplices){													//	for j=0 to m-1 {
		if(tArray[curIndex].marked && tArray[curIndex].simplex.size() == 0){							//		if sj is marked and T[j] is empty{
			intervals[simp.first.size()].push_back(std::make_pair(kSimplices[curIndex].second, maxEpsilon));		//			k = dim sj; Lk = Lk U {(deg sj, inf)}
		}
	}	*/ 
	std::cout << "return" << std::endl;
	return ret;
}


std::set<unsigned> persistencePairs::removePivotRows(int currentIndex, std::pair<std::set<unsigned>, double> simp, std::vector<std::vector<std::vector<unsigned>>> dimBoundaries){
	//def RemovePivotRows(s):
	//k = dim s; d = delk*s
	//Remove unmarked items in d;
	//while (d!= 0){
	//	i = maxindex d;
	//	if T[i] is empty, break;
	//	Let q be the coefficient of si in T[i]
	//	d = d-q^-1*T[i]
	//}
	//return d;
	
	//In basic terms: 
	//	Look at the boundary matrix and determine if the current indices are a pivot
	//	This should be the columns of dimBoundaries[k]
	//		if this is a pivot, return a null set; this will create a k-1 betti
	//			for the lowest index k-1 simplex (i.e. [a,b,c] -> [a,b]; [a,b] -> a;
	//			(We can just cut off the last index in the simplex set)
	//		if this is not a pivot, we need to mark the row (as it will be used later)
	//			and return nothing :)
	
	
	//Determine if the current index is a pivot
	//if(
	
	
	
	int k = simp.first.size();
	std::vector<std::vector<unsigned>> d = dimBoundaries[k];
	
	
	
	while(d.size() > 0){
		int i = *simp.first.rbegin();
		//if(
		
	}
	
	
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
pipePacket persistencePairs::runPipe(pipePacket inData){

	std::vector<int> bettiNumbers (dim,0);
	std::vector<float> lifeSpans[dim];
	
	std::vector<std::vector<std::pair<double, double>>> intervals;
	
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges = inData.workData.complex->getAllEdges(maxEpsilon);
	
	//Retrieve
	auto local_weights = inData.weights;
	
	double epsilon = 0;
	
	//Compute the persistence pairs as detailed by Zomorodian
	
	computeIntervals(edges);
		
	return inData;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void persistencePairs::outputData(pipePacket inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_bettis_output.csv");
	
	file << inData.bettiOutput;
	
	file.close();
	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool persistencePairs::configPipe(std::map<std::string, std::string> configMap){
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

