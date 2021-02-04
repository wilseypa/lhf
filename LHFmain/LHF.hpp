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

	int nprocs = 1, id = 0;

	void testFunc(int num1);
	void myprint(void);
	void outputBettis(std::map<std::string, std::string>, pipePacket &);
	void runPipeline(std::map<std::string, std::string>, pipePacket &);
	void processDataWrapper(std::map<std::string, std::string>, pipePacket &);
	std::vector<bettiBoundaryTableEntry> processParallel(std::map<std::string, std::string>, std::vector<unsigned>&, std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>>&, int = 0);
	std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string>, pipePacket &, bool = true);
	std::vector<bettiBoundaryTableEntry> processUpscaleWrapper(std::map<std::string, std::string>, pipePacket &);
};
extern "C" {

	LHF void outputBettis(std::map<std::string, std::string>, pipePacket &);
	LHF void runPipeline(std::map<std::string, std::string>, pipePacket &);
	LHF void processDataWrapper(std::map<std::string, std::string>, pipePacket &);
	LHF std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string> args, pipePacket &wd){return processIterUpscale(args, wd, true);}
	LHF std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string>, pipePacket &, bool);
	LHF std::vector<bettiBoundaryTableEntry> processUpscaleWrapper(std::map<std::string, std::string>, pipePacket &);
	
	//wrapper for c++ class

	LHF void wrapper_init_class(int id);
	LHF void wrapper_free_class();
	
	//wrapper versions?

	// LHF void wrapper_outputBettis(std::map<std::string, std::string>, pipePacket &);
	// LHF void wrapper_runPipeline(std::map<std::string, std::string>, pipePacket &);
	// LHF void wrapper_processDataWrapper(std::map<std::string, std::string>, pipePacket &);
	// LHF std::vector<bettiBoundaryTableEntry> wrapper_processIterUpscale(std::map<std::string, std::string> args, pipePacket &wd){return processIterUpscale(args, wd, true);}
	// LHF std::vector<bettiBoundaryTableEntry> wrapper_processIterUpscale(std::map<std::string, std::string>, pipePacket &, bool);
	// LHF std::vector<bettiBoundaryTableEntry> wrapper_processUpscaleWrapper(std::map<std::string, std::string>, pipePacket &); 
}