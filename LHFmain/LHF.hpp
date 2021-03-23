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
		return lhs.bettiDim < rhs.bettiDim || (lhs.bettiDim == rhs.bettiDim && lhs.birth > rhs.birth);
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
	std::vector<bettiBoundaryTableEntry> processParallelWrapper(std::map<std::string, std::string>, pipePacket &, bool = true);
	std::vector<bettiBoundaryTableEntry> processDistributedWrapper(std::map<std::string, std::string>, pipePacket &);

};



extern "C" {
	//Handle Betti Return Structure allocation	
	typedef struct bRetStructure {
		int dim;
		double birth,death;
	} BRET;
	
	
	void free_bRet(BRET *b){
		free(b);
	}
	
	typedef struct bWrapStructure{
		int size;
		BRET* ret;
	} BRAP;
	
	
	void testFunc(int num1, char* st) { std::cout << "Test: " << num1 << std::endl; std::cout << "\t" << st << std::endl;};
	void pyRunWrapper(const int, char*, const double *);
	BRAP* pyRunWrapper2(const int, char*, const double *);
}
