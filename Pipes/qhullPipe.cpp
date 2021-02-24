
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
#include "qhullPipe.hpp"
#include "utils.hpp"

// basePipe constructor
qhullPipe::qhullPipe(){
	pipeType = "qhullPipe";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
void qhullPipe::runPipe(pipePacket &inData){
	
    RboxPoints eg("100");
    Qhull q(eg, "");
    QhullFacetList facets= q.facetList();
    std::cout << facets;

	ut.writeDebug("qhullPipe", "\tSuccessfully Executed pipe");
	return;
}


// configPipe -> configure the function settings of this pipeline segment
bool qhullPipe::configPipe(std::map<std::string, std::string> &configMap){
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

	configured = true;
	ut.writeDebug("qhullPipe","Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
void qhullPipe::outputData(pipePacket &inData){
	std::ofstream file;
	file.open("output/" + pipeType + "_output.csv");
	
	
	
	file.close();
	return;
}

