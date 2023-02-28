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



template<class nodeType>
class LHF {
  private:
	
  public:
	
	int nprocs = 1, id = 0;

	void testFunc(int num1) { std::cout << num1 << std::endl;};
	void myprint(void);
	void outputBettis(std::map<std::string, std::string>, pipePacket<nodeType> &);
	void runPipeline(std::map<std::string, std::string>, pipePacket<nodeType> &);
	void runPreprocessor(std::map<std::string, std::string> &, pipePacket<nodeType> &);
	std::vector<bettiBoundaryTableEntry> processParallel(std::map<std::string, std::string>, std::vector<unsigned>&, std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>>&, std::vector<std::vector<double>>&, int = 0);
	std::vector<bettiBoundaryTableEntry> processParallelWrapper(std::map<std::string, std::string>, pipePacket<nodeType> &, bool = true);
	std::vector<bettiBoundaryTableEntry> processDistributedWrapper(std::map<std::string, std::string>, pipePacket<nodeType> &);

};

//Explicit Template Class Instantiation
template class LHF<simplexNode>;
template class LHF<alphaNode>;
template class LHF<witnessNode>;

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

	//### pipepacket ###//

	// void free_pRet(PRET *b){
	// 	free(b);
	// }

	typedef struct pipeWrapStructure{
		int size_betti; //multiple sizes?
		int LHF_size;
		int LHF_dim;
		int workData_size;
		BRET* BettiTable;
		double* inputData;
		double* distMatrix;
		double* workData;
		unsigned* centroidLabels;
		char* stats;
		char* runLog;
		char* ident;
	} PRAP;

	//##################//
	void testFunc(int num1, char* st) { std::cout << "Test: " << num1 << std::endl; std::cout << "\t" << st << std::endl;};
	void pyRunWrapper(const int, char*, const double *);
	PRAP* pyRunWrapper2(int, char *, const double *);
}
