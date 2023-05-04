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
    pipePacket<nodeType> *processpyLHFWrapper(std::map<std::string, std::string> &, std::vector<std::vector<double>> &); 

};

//Explicit Template Class Instantiation
template class LHF<simplexNode>;
template class LHF<alphaNode>;
template class LHF<witnessNode>;

extern "C" {
    
    /******************************** NEW INTERFACE *****************************/
    typedef struct retBettiTable{
	/**
     * @brief BettiTable return structure for the python interface; encodes resultant PIs and generators.
     * 
     * This structure captures the serialized betti boundary table for passing to the pyLHF wrapper. These
     *      structures must match the ctypes definitions in LHF.py
     * 
     */
        int dim;                        // Homology class of the PI
        double birth, death;            // Birth and death times of the PI
        int boundarySize;               // # of entries (generators) for the boundary
        unsigned* boundaryEntries;      // Generators of the boundary (set)
        
    } bettiTableWrap;
    
    
    typedef struct retPipePacket{
	/**
     * @brief Return structure for the python interface; encodes the entire pipePacket object
     * 
     * This structure captures the serialized pipePacket with size information for passing
     *      to the pyLHF wrapper. These structures must match the ctypes definitions in LHF.py.
     * 
     */
    	int size_betti;                             // # of entries (vectors) in the betti table
		int LHF_size;                               // # of entries (vectors) in the LHF data set
		int LHF_dim;                                // Dimension of entries (vectors) in the LHF data set
		int workData_size;                          // # of entries (vectors) in the working data set
		bettiTableWrap* bettiTable  = nullptr;      // Pointer to the betti table structure
		double* inputData           = nullptr;      // Pointer to serialized input data
		double* distMatrix          = nullptr;      // Pointer to distance matrix
		double* workData            = nullptr;      // Pointer to working data
		unsigned* centroidLabels    = nullptr;      // Pointer to centroidLabels
		const char* stats;                          // Pointer to string stats
		const char* runLog;                         // Pointer to string runLog
		const char* ident;                          // Pointer to string identity
    } pipeWrap;
        
    void free_bettiWrap(bettiTableWrap *b){
	/**
     * @brief Free the bettiTable structure after enumerating the data in pyLHF.
     */
        free(b); 
        return; 
    }
    
    void free_pipeWrap(pipeWrap *p){ 
	/**
     * @brief Free the pipepacket structure after enumerating the data in pyLHF; if the bettiTable
     *      has not been previously freed this function will also call free_bettiWrap.
     */
        if(p->size_betti > 0){
            free_bettiWrap(p->bettiTable);
        }
        delete(p); 
        return; 
    }
    
    pipeWrap* pyLHFWrapper(int, char *, const double *);
       
    
}
