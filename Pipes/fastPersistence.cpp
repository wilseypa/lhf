/*
 * persistencePairs hpp + cpp extend the basePipe class for calculating the 
 * persistence pairs numbers from a complex
 * 
 */

#include <string>
#include <chrono>
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
#include "fastPersistence.hpp"


// basePipe constructor
fastPersistence::fastPersistence(){
	pipeType = "FastPersistence";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//	FastPersistence: For computing the persistence pairs from simplicial complex:
//		1. See Bauer-19 for algorithm/description
pipePacket fastPersistence::runPipe(pipePacket inData){
	
	std::string bettis = "";
	
	//Get all edges for the simplexArrayList or simplexTree
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> edges = inData.complex->getAllEdges(maxEpsilon);
	
	if(edges.size() < 2){
		ut.writeLog("FastPersistence","Failed to provide the fastPersistence pipe with edges from complex");
		return inData;
	}
	
	std::vector<std::pair<double,double>> temp;
	std::vector<std::vector<std::pair<double,double>>> ret;
	for(int i = 0; i < dim; i++){
		ret.push_back(temp);
	}
	
	//Some notes on fast persistence:
	
	//	-Vectors need to be stored in a lexicograhically ordered set of decreasing (d+1)-tuples (e.g. {3, 1, 0})
	//		-These vectors are replaced with their indices (e.g. {{2,1,0} = 0, {3,1,0} = 1, {3,2,0} = 2, etc.})
	
	//	-Reduction matrix (V) stores collection of non-zero entries for each column (indexed)
	
	
	//Start a timer for physical time passed during the pipe's function
	auto startTime = std::chrono::high_resolution_clock::now();
	
	//Get all dim 0 persistence intervals
		//Kruskal's minimum spanning tree algorithm
	
	std::vector<std::pair<std::set<unsigned>,double>> mst;
	std::set<unsigned> conSet;
	std::set<unsigned> wset;
		
	for(auto edge : edges[1]){
		std::cout << edge.second << "\t";
		if((wset = ut.setIntersect(edge.first, conSet, false)).size() < 2){
			mst.push_back(edge);
			for(auto i = edge.first.begin(); i != edge.first.end(); i++){
				conSet.insert(*i);
			}
			
		}
	}
	
	std::cout << "D0: " << std::endl;
	for(auto z : mst){
		std::cout << z.second << "\t";
		ut.print1DVector(z.first);
	}
	std::cout << std::endl;
	
	//For higher dimensional persistence intervals
	//	
	
	
	//
	
	
	//Stop the timer for time passed during the pipe's function
	auto endTime = std::chrono::high_resolution_clock::now();
	
	//Calculate the duration (physical time) for the pipe's function
	std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
	
	//Output the time and memory used for this pipeline segment
	ut.writeDebug("persistence","Bettis executed in " + std::to_string(elapsed.count()/1000.0) + " seconds (physical time)");;
	
	//Print the bettis
	if(debug)
		std::cout << std::endl << bettis << std::endl;
		
	inData.bettiOutput = bettis;
		
	return inData;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void fastPersistence::outputData(pipePacket inData){
	std::ofstream file;
	if(fnmod.size() > 0)
		file.open("output/"+pipeType+"_bettis_output"+fnmod+".csv");
	else
		file.open("output/" + pipeType + "_bettis_output.csv");
	
	file << inData.bettiOutput;
	
	file.close();
	
	file.open("output/tArray.csv");
	
	file << "ti,Dim,Marked,Birth,Death,Simplex\n";
	for(auto tStruct : tArray){
		file << tStruct.ti << "," << tStruct.simplex.size()-1 << "," << tStruct.marked << "," << tStruct.birth << "," << tStruct.death << ",";
		for(auto index : tStruct.simplex)
			file << index << " ";
		file << "\n";
	}
	file.close();
	
	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool fastPersistence::configPipe(std::map<std::string, std::string> configMap){
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
	
	pipe = configMap.find("twist");
	if(pipe != configMap.end())
		twist = configMap["twist"];
	else return false;
	
	pipe = configMap.find("fn");
	if(pipe != configMap.end())
		fnmod = configMap["fn"];
	
	pipe = configMap.find("complexType");
	if(pipe != configMap.end() && configMap["complexType"] == "indSimplexTree")
		alterPipe = true;
		
	configured = true;
	ut.writeDebug("fastPersistence","Configured with parameters { dim: " + configMap["dimensions"] + ", twist: " + twist + ", complexType: " + configMap["complexType"] + ", eps: " + configMap["epsilon"]);
	ut.writeDebug("fastPersistence","\t\t\t\tdebug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

