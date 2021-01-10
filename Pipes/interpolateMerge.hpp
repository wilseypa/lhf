#pragma once

// Header file for buldingValidSimplicialComplexes class - see buildComplexes.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class interpolateMergePipe : public basePipe {
	private:
 	public:
		int dim;
	    interpolateMergePipe();
	    void runPipe(pipePacket &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
		void outputData(pipePacket&);
};
