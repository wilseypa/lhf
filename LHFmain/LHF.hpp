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
    void* processpyLHFWrapper(std::map<std::string, std::string> &, std::vector<std::vector<double>> &); 

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
    	int size_betti;                 // # of entries (vectors) in the betti table
		int LHF_size;                   // # of entries (vectors) in the LHF data set
		int LHF_dim;                    // Dimension of entries (vectors) in the LHF data set
		int workData_size;              // # of entries (vectors) in the working data set
		bettiTableWrap* bettiTable;     // Pointer to the betti table structure
		double* inputData;              // Pointer to serialized input data
		double* distMatrix;             // Pointer to distance matrix
		double* workData;               // Pointer to working data
		unsigned* centroidLabels;       // Pointer to centroidLabels
		const char* stats;                    // Pointer to string stats
		const char* runLog;                   // Pointer to string runLog
		const char* ident;                    // Pointer to string identity
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
        if(p->bettiTable != nullptr){
            free_bettiWrap(p->bettiTable);
        }
        free(p); 
        return; 
    }
    
    pipeWrap* pyLHFWrapper(int, char *, const double *);
        
        
        
        
    /******************************** OLD *****************************/
        
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
