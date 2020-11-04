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

/** LHF Library Class
 * 
 *  Implements shared library functions for general use. These interfaces are the general approaches to
 *      using several LHF modes including PPH, streaming, and fastPersistence.
 */
class LHF {
  private:
	
  public:
	
	/** @brief Total number of MPI processes for execution
	 *  
	 *  Stores the total number of MPI processes executing on the distributed set of partitions.
	 */
	int nprocs;

	/** @brief ID of the current MPI process
	 * 
	 *  Stores the current MPI process ID for executing on the distributed partitions.
	 */
	int id;


	/** outputBettis LHFLib function
	 * 
	 *  @brief Output pipePacket persistence intervals to file
	 * 
	 *  OutputBettis writes persistence intervals contained in the passed pipepacket
	 *      into argument driven file ("outputFile")
	 * 
	 *  @param args std::map<std::string,std::string> of arguments
	 *  @param &wD pipePacket with persistence intervals after execution
	 */
	void outputBettis(std::map<std::string, std::string>, pipePacket &);
	
	
	/** runPipeline LHFLib function
	 * 
	 *  @brief Run the primary pipeline components
	 * 
	 *  Run the primary pipeline components set through ["pipeline"] argument mapping. Does not
	 *      include the preprocessor function in this function. Changes are stored
	 *      into referenced pipePacket object.
	 * 
	 *  @param args std::map<std::string,std::string> of arguments
	 *  @param &wD pipePacket after data has been read (unless using a streaming type)
	 */
	void runPipeline(std::map<std::string, std::string>, pipePacket &);
	
	
	/** processDataWrapper LHFLib function
	 * 
	 *  @brief Run the preprocessor pipeline function
	 * 
	 *  Run the preprocessor for a configured function using the ["preprocessor"] argument mapping. Changes are stored
	 *     into referenced pipePacket object.
	 * 
	 *  @param args std::map<std::string,std::string> of arguments
	 *  @param &wD pipePacket after data has been read (unless using a streaming type)
	 */
	void processDataWrapper(std::map<std::string, std::string>, pipePacket &);
	
	
	/** processIterUpscale LHFLib function
	 * 
	 *  @brief Run partitioned persistent homology
	 * 
	 *  Main entrypoint for performing partitioned persistent homology with upscaling and/or iterative reduction of
	 *     input point cloud. Changes are stored into referenced pipePacket object.
	 * 
	 *  @param args std::map<std::string,std::string> of arguments
	 *  @param &wD pipePacket after data has been read (unless using a streaming type)
	 */
	std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string> args, pipePacket &wd){return processIterUpscale(args, wd, true);}
	
	
	/** processIterUpscale LHFLib function
	 * 
	 *  @brief Run partitioned persistent homology with option to repartition
	 * 
	 *  Main entrypoint for performing partitioned persistent homology with upscaling and/or iterative reduction of
	 *     input point cloud. Changes are stored into referenced pipePacket object.
	 * 
	 *  @param args std::map<std::string,std::string> of arguments
	 *  @param &wD pipePacket after data has been read (unless using a streaming type)
	 *  @param runPartition boolean to determine if partitioning of the input data is required
	 */
	std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string>, pipePacket &, bool);
	
	
	/** processUpscaleWrapper LHFLib function
	 * 
	 *  @brief Run partitioned persistent homology using MPI for partition distribution
	 * 
	 *  Main entrypoint for performing partitioned persistent homology with upscaling and/or iterative reduction of
	 *     input point cloud. Partitions are distributed to specified nodes with OpenMPI to provide distributed speedup.
	 *     Changes are stored into referenced pipePacket object.
	 * 
	 *  @param args std::map<std::string,std::string> of arguments
	 *  @param &wD pipePacket after data has been read (unless using a streaming type)
	 */
	std::vector<bettiBoundaryTableEntry> processUpscaleWrapper(std::map<std::string, std::string>, pipePacket &);

};
