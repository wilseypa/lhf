 /*
 *Build Complexes for all the valid_d_simplexmeshes
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
#include <set>
#include <algorithm>
#include "distBuildSCExecutePH.hpp"
#include "utils.hpp"

// basePipe constructor
distBuildSCExecutePHPipe::distBuildSCExecutePHPipe(){
	pipeType = "distBuildSCExecutePH";
	return;
}

std::vector<bettiBoundaryTableEntry>  vscPersistence::lightPersistence(){

	std::vector<bettiBoundaryTableEntry> bettiTable;

	return bettiTable;
}	// Complex 
// runPipe -> Run the configured functions of this pipeline segment
void distBuildSCExecutePHPipe::runPipe(pipePacket &inData){
int counter =0;
std::vector<std::set<simplexNode_P, cmpByWeight>> simplicialComplex;          // Complex
for(auto validSC : inData.VUdsimplexes)
{   
	simplicialComplex.clear();          // Complex

	simplicialComplex =  inData.complex->buildValidSimplicialComplex(validSC,inData.inputData.size());
	if(counter%1000==0)
		std::cout<<counter;
	counter++;
	double aEW=0;
        for(auto edge : simplicialComplex[1])
		aEW = aEW + (*edge).weight;
	aEW/= (double)simplicialComplex[1].size();

	auto vPH = vscPersistence(simplicialComplex);
	std::vector<bettiBoundaryTableEntry> bettiTable = vPH.lightPersistence();
	

	bettiTables BT =  bettiTables(bettiTable,(*(simplicialComplex[1].rbegin()))->weight,aEW) ;

        inData.bTbs.insert(BT);
}	
return;
	
}



// configPipe -> configure the function settings of this pipeline segment
bool distBuildSCExecutePHPipe::configPipe(std::map<std::string, std::string> &configMap){
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
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		epsilon = std::atof(configMap["epsilon"].c_str());
	else return false;
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	else return false;
	
        configured = true;	
	ut.writeDebug("distBuildSCExecutePHPipe","Configured with parameters { dim: " + std::to_string(dim) + " , debug: " + strDebug + ", outputFile: " + outputFile);
	
	return true;
}



void distBuildSCExecutePHPipe::outputData(pipePacket &inData){
	return;	
}	
