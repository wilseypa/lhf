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

	void testFunc(int num1) { std::cout << num1 << std::endl;};
	void myprint(void);
	void outputBettis(std::map<std::string, std::string>, pipePacket &);
	void runPipeline(std::map<std::string, std::string>, pipePacket &);
	void processDataWrapper(std::map<std::string, std::string>, pipePacket &);
	std::vector<bettiBoundaryTableEntry> processParallel(std::map<std::string, std::string>, std::vector<unsigned>&, std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>>&, int = 0);
	std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string>, pipePacket &, bool = true);
	std::vector<bettiBoundaryTableEntry> processUpscaleWrapper(std::map<std::string, std::string>, pipePacket &);
	
};



extern "C" {
	void testFunc(int num1) { std::cout << "Test: " << num1 << std::endl;};
	void pyRunWrapper(std::map<std::string, std::string>, std::vector<std::vector<double>>);
}
