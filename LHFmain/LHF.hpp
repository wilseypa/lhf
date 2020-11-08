#pragma once

#include "readInput.hpp"
#include "argParser.hpp"
#include "basePipe.hpp"
#include "writeOutput.hpp"
#include "pipePacket.hpp"
#include "preprocessor.hpp"
#include "utils.hpp"

struct sortBettis{
	bool operator()(bettiBoundaryTableEntry lhs, bettiBoundaryTableEntry rhs){
		return lhs.bettiDim < rhs.bettiDim || (lhs.bettiDim == rhs.bettiDim && lhs.death < rhs.death);
	}
};

class LHF {
  private:
	
  public:
	int nprocs,id;

	void outputBettis(std::map<std::string, std::string>, pipePacket &);
	void runPipeline(std::map<std::string, std::string>, pipePacket &);
	void processDataWrapper(std::map<std::string, std::string>, pipePacket &);
	std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string>, pipePacket &, bool = true);
	std::vector<bettiBoundaryTableEntry> processUpscaleWrapper(std::map<std::string, std::string>, pipePacket &);
};
