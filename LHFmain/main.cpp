#include "mpi.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <thread>
#include "readInput.hpp"
#include "argParser.hpp"
#include "basePipe.hpp"
#include "writeOutput.hpp"
#include "pipePacket.hpp"
#include "preprocessor.hpp"
#include "utils.hpp"

using namespace std;


void runPipeline(std::map<std::string, std::string> args, pipePacket* wD){
	auto *ws = new writeOutput();

	// Begin processing parts of the pipeline
	// DataInput -> A -> B -> ... -> DataOutput
	// Parsed by "." -> i.e. A.B.C.D
	auto pipe = args.find("pipeline");
	if(pipe != args.end()){
		auto pipeFuncts = std::string(args["pipeline"]);
		auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;
		
		//For each '.' separated pipeline function (count of '.' + 1 -> lim)
		for(unsigned i = 0; i < lim; i++){
			auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
			pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);
			
			//Build the pipe component, configure and run
			auto *bp = new basePipe();
			auto *cp = bp->newPipe(curFunct, args["complexType"]);
			
			//Check if the pipe was created and configure
			if(cp != 0 && cp->configPipe(args)){
				//Run the pipe function (wrapper)
				*wD = cp->runPipeWrapper(*wD);
			} else {
				cout << "LHF : Failed to configure pipeline: " << args["pipeline"] << endl;
			}
		}
	}
	//If the pipeline was undefined...
	else {
		cout << "LHF : Failed to find a suitable pipeline, exiting..." << endl;
		return;
	}
	
	//Output the data using writeOutput library
	pipe = args.find("outputFile");
	if(pipe != args.end()){
		if (args["outputFile"] == "console"){
			//ws->writeConsole(wD);
		} else {
			ws->writeStats(wD->stats, args["outputFile"]);
			ws->writeBarcodes(wD->bettiOutput, args["outputFile"]);
			
		}
	}
	
	return;
}



void processDataWrapper(std::map<std::string, std::string> args, pipePacket* wD){
    
	//Start with the preprocessing function, if enabled
	auto pre = args["preprocessor"];
	if(pre != ""){
		auto *preprocess = new preprocessor();
		auto *prePipe = preprocess->newPreprocessor(pre);
		
		if(prePipe != 0 && prePipe->configPreprocessor(args)){
			*wD = prePipe->runPreprocessorWrapper(*wD);
		} else {
			cout << "LHF : Failed to configure pipeline: " << args["pipeline"] << endl;
		}
	}
	
	runPipeline(args, wD);
		
	return;
}	

void processReducedWrapper(std::map<std::string, std::string> args, pipePacket* wD){
    
    std::vector<bettiBoundaryTableEntry> mergedBettiTable;
    
	//Start with the preprocessing function, if enabled
	auto pre = args["preprocessor"];
	if(pre != ""){
		auto *preprocess = new preprocessor();
		auto *prePipe = preprocess->newPreprocessor(pre);
		
		if(prePipe != 0 && prePipe->configPreprocessor(args)){
			*wD = prePipe->runPreprocessorWrapper(*wD);
		} else {
			cout << "LHF : Failed to configure pipeline: " << args["pipeline"] << endl;
		}
	}
	
	utils ut;
	
	auto maxRadius = ut.computeMaxRadius(std::atoi(args["clusters"].c_str()), wD->originalData, wD->fullData, wD->originalLabels);
	auto centroids = wD->originalData;
	std::cout << "Using maxRadius: " << maxRadius << std::endl;
	
	auto partitionedData = ut.separatePartitions(2*maxRadius, wD->originalData, wD->fullData, wD->originalLabels);
	
	std::cout << "Partitions: " << partitionedData.second.size() << std::endl << "Counts: ";
	for(unsigned z = 0; z < partitionedData.second.size(); z++){
		if(partitionedData.second[z].size() > 0){
			std::cout << "Running Pipeline with : " << partitionedData.second[z].size() << " vectors" << std::endl;
			wD->originalData = partitionedData.second[z];
			
			
			runPipeline(args, wD);
			
			wD->complex->clear();
			
			//Map partitions back to original point indexing
			ut.mapPartitionIndexing(partitionedData.first[z], wD->bettiTable);
			
			for(auto betEntry : wD->bettiTable)
				mergedBettiTable.push_back(betEntry);
			
		} else 
			std::cout << "skipping" << std::endl;
	}
	
	std::cout << "Full Data: " << centroids.size() << std::endl;
	if(centroids.size() > 0){
		std::cout << "Running Pipeline with : " << centroids.size() << " vectors" << std::endl;
		wD->originalData = centroids;
		runPipeline(args, wD);
		
		wD->complex->clear();
	} else 
		std::cout << "skipping" << std::endl;
	
	std::cout << std::endl << "_______BETTIS_______" << std::endl;
	
	for(auto a : wD->bettiTable){
		std::cout << a.bettiDim << ",\t" << a.birth << ",\t" << a.death << ",\t";
		ut.print1DVector(a.boundaryPoints);
	}
	
		
	return;
}	





void processUpscaleWrapper(std::map<std::string, std::string> args, pipePacket* wD){
    
	//Start with the preprocessing function, if enabled
	auto pre = args["preprocessor"];
	if(pre != ""){
		auto *preprocess = new preprocessor();
		auto *prePipe = preprocess->newPreprocessor(pre);
		
		if(prePipe != 0 && prePipe->configPreprocessor(args)){
			*wD = prePipe->runPreprocessorWrapper(*wD);
		} else {
			cout << "LHF : Failed to configure pipeline: " << args["pipeline"] << endl;
		}
	}
	
	
	//Notes from Nick:
	//	After preprocessor, need to separate the data into bins by index (partition)
	//		-Look at each point, label retrieved from preprocessor
	//			-bin points and send to a respective MPI node / slave
	utils ut;
	auto partitionedData = ut.separatePartitions(std::atoi(args["clusters"].c_str()), wD->fullData, wD->originalLabels);
	
	//	Each node/slave will process at least 1 partition
	//		NOTE: the partition may contain points outside that are within 2*Rmax
	std::cout << "Partitions: " << partitionedData.second.size() << std::endl << "Counts: ";
	
	for(auto a : partitionedData.second){
		std::cout << a.size() << "\t";
	}
	std::cout << std::endl;
	
	// Once we have X data sets to process, along with the original centroid dataset
	//		Distribute work to slaves (run fastPersistence on given data)
	//		Retrieve results
	//		Merge Barcodes	
	
	
	
	
	
	
	
	
	
	
	
	// OLD----- (but keeep for reference later)
	
	/*do{
		if(wD->boundaries.size() > 0){
			
			std::vector<std::thread> threads;
			
			//Spin off a pipeline for each boundary...
			for(auto bound : wD->boundaries){
				std::thread curThread([&]{runPipeline(args,wD);});
				//curThread.detach();
				threads.push_back(std::move(curThread));
			}
			for(int a = 0; a < threads.size(); a++){
				threads[a].join();
			}
		} else {
			runPipeline(args, wD);
		}
	} while (wD->boundaries.size() > 0);
		
	return;*/
}	

int main(int argc, char* argv[]){
	//	Steps to compute PH with preprocessing:
	//
	//		1.  Preprocess the input data
	//			a. Store the original dataset, reduced dataset, cluster indices
	//
	//		2.	Build a 2-D neighborhood graph
	//			a. Weight each edge less than epsilon
	//			b. Order edges (min to max)
	//
	//		3.  Rips filtration
	//			a. Build higher-level simplices (d > 2)
	//				- Retain the largest of edges as weight (when simplex forms)
	//			b. Build list of **relevant epsilon values
	//
	//		4.	Betti calculation
	//			a. For each relevant epsilon value
	//				i.	Create Boundary Matrix
	//					- if pn[i,j] < n, set to 1
	//					- else set to 0
	//				ii.	Compute RREF
	//				iii.Store constituent boundary points
	//
	//		5.	Upscaling
	//			a. Upscale around boundary points
	//			b. For each upscaled boundary, repeat steps 2-4
	//
	//
	//	**relevant refers to edge weights, i.e. where a simplex of any
	//			dimension is created (or merged into another simplex)
	
	
	//Define external classes used for reading input, parsing arguments, writing output
	auto *rs = new readInput();
    auto *ap = new argParser();
    
    //Parse the command-line arguments
    auto args = ap->parse(argc, argv);
    
    //Determine what pipe we will be running
    ap->setPipeline(args);
    
    for(auto z : args)
		std::cout << z.first << "\t" << z.second << std::endl;
    
	//Create a pipePacket (datatype) to store the complex and pass between engines
    auto *wD = new pipePacket(args, args["complexType"]);	//wD (workingData)
	
	if(args["pipeline"] != "slidingwindow"){
		//Read data from inputFile CSV
		wD->originalData = rs->readCSV(args["inputFile"]);
		wD->fullData = wD->originalData;
	}
	//If data was found in the inputFile
	if(wD->originalData.size() > 0 || args["pipeline"] == "slidingwindow"){
		
		//Add data to our pipePacket
		wD->originalData = wD->originalData;
		
		if(args["upscale"] == "true"){
			processUpscaleWrapper(args, wD);
		} else if (args["mode"] == "reduced"){
			processReducedWrapper(args,wD);			
		} else {
			processDataWrapper(args, wD);
		}
	} else {
		ap->printUsage();
	}
	
	
	
    return 0;
}
