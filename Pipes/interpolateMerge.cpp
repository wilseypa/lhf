/*
 *Interpolate and merge Betties for all diffrent executions
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
#include "interpolateMerge.hpp"
#include "utils.hpp"

// basePipe constructor
interpolateMergePipe::interpolateMergePipe(){
	pipeType = "interpolateMerge";
	return;
}


void interpolateMergePipe::runPipe(pipePacket &inData){
//std::vector<bettiBoundaryTableEntry> finalbettiTable;
std::ofstream file("bettisequence.csv");
/*for(auto entry : inData.bTbs){
 //   file<<entry.bettiTab.size()<<","<<entry.numpts<<","<<entry.maxEpsilon;
 //   file<<endl;
    for(auto a : entry.bettiTab){
        file<<a.bettiDim<<","<<a.birth<<","<<a.death<<std::endl;
        bool found = false;
        for(auto x : finalbettiTable)
        {
           if(a.bettiDim==x.bettiDim&&a.birth==x.birth&&a.death==x.death&&a.boundaryPoints==x.boundaryPoints)
               {
                found = true;
                break;
               }
        }
        if(!found)
            finalbettiTable.push_back(a);
    }
 //   file<<endl;
}
*/
file<<std::endl<<"Final Betti Table"<<std::endl;
for(auto a : inData.bettiTable){
        file<<a.bettiDim<<","<<a.birth<<","<<a.death;
        for(auto k : a.boundaryPoints)
            file<<k<<",";
    file<<std::endl;
}
file.close();



}

// configPipe -> configure the function settings of this pipeline segment
bool interpolateMergePipe::configPipe(std::map<std::string, std::string> &configMap){
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
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
        configured = true;	
	ut.writeDebug("interpolateMergePipe","Configured with parameters { dim: " + std::to_string(dim) + " , debug: " + strDebug + ", outputFile: " + outputFile);
	
	return true;
}


void interpolateMergePipe::outputData(pipePacket &inData){
	return;		
}
