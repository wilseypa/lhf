#pragma once

// Header file for enumerateUniqueValidSimplicialComplexes class - see enumerateUVSC.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class enumerateUVSCPipe : public basePipe {
	
 	public:
	    int dim;
	    enumerateUVSCPipe();
	    void runPipe(pipePacket &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
            void outputData(pipePacket&);
};

