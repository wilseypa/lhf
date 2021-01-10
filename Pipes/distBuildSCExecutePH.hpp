#pragma once

// Header file for buldingValidSimplicialComplexes class - see buildComplexes.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "argParser.hpp"

class vscPersistence{

	private:
        	std::vector<std::set<simplexNode_P, cmpByWeight>> simplicialComplex;          // Complex 
	public:
	vscPersistence();
	vscPersistence(std::vector<std::set<simplexNode_P, cmpByWeight>> sc) {
		simplicialComplex = sc;
	};		// Complex 
	std::vector<bettiBoundaryTableEntry>  lightPersistence();
};
class distBuildSCExecutePHPipe : public basePipe {
	private:
 	public:
           int dim;
	   distBuildSCExecutePHPipe();
           void runPipe(pipePacket &inData);
	   bool configPipe(std::map<std::string, std::string> &configMap);
           void outputData(pipePacket&);
};
