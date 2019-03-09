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

void persistencePairs::computeIntervals(std::vector<std::pair<std::set<unsigned>, double>> kSimplices, std::vector<std::vector<std::pair<double,double>>> &intervals){
	struct tArrayEntry_t{
		bool marked = false;
		int maxIndex;
		std::set<unsigned> simplex;
	};
		
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
	
	//Create our T-array
	tArrayEntry_t tArray[kSimplices.size()];
	
	std::cout << "total simplices: " << kSimplices.size() << std::endl;
	
	
	int curIndex = 0; //Track our index in the tArray
	
	//for k=0 to dim(K) Lk = 0;
	for(auto simp : kSimplices){													//	for j=0 to m-1 {
		std::set<unsigned> d = removePivotRows(simp);								//		d = RemovePivotrows(sj)
		if(d.size() == 0)
			tArray[curIndex].marked = true;											//		if(d==0) Mark sj
		else {
			unsigned i = *d.rbegin();												//			i = maxindexd; k = dim si;
			auto k = kSimplices[i].first.size();									// !! need to get kSimplices[i].size in a diff manner
			tArray[i].maxIndex = curIndex;
			tArray[i].simplex = d;													//			store j and d in T[i];
			
			intervals[simp.first.size()].push_back(std::make_pair(kSimplices[i].second, kSimplices[curIndex].second));
			
																					//			Lk = Lk U {(deg si, deg sj)}
																					//NOTE: Degree is the weight of the simplex
		}
		curIndex++;
	}
	
	curIndex = 0;
	for(auto simp : kSimplices){													//	for j=0 to m-1 {
		if(tArray[curIndex].marked && tArray[curIndex].simplex.size() == 0){							//		if sj is marked and T[j] is empty{
			intervals[simp.first.size()].push_back(std::make_pair(kSimplices[curIndex].second, maxEpsilon));		//			k = dim sj; Lk = Lk U {(deg sj, inf)}
		}
	}	
	
	return;
}


std::set<unsigned> persistencePairs::removePivotRows(std::pair<std::set<unsigned>, double> &simp){
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
	std::vector<std::pair<std::set<unsigned>,double>> flattenedEdges;
	
	for(auto a : edges){
		for(auto z : a){
			flattenedEdges.push_back(z);
		}
	}
	
	
	//Retrieve
	auto local_weights = inData.weights;

	std::string bettiOutput[dim+1];
	std::fill_n(bettiOutput, dim+1, "Epsilon\t\tDim\tBD\tBetti\trank\tnullity\tlRank\tlNul\n");
	
	double epsilon = 0;
	
	//Compute the persistence pairs as detailed by Zomorodian
	
	computeIntervals(flattenedEdges, intervals);
		
	if(debug){
		for(int d = 0; d < dim; d++){
			std::cout << "______Intervals d = " << d << "______" << std::endl;
			std::cout << "dim\tbirth\tdeath" << std::endl;
			for (auto a : intervals[d])
				std::cout << d << "\t" << a.first << "\t" << a.second << std::endl;
			std::cout << std::endl;
		}
	}
	
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

