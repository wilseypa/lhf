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
			auto *cp = bp->newPipe(curFunct);
			
			//Check if the pipe was created and configure
			if(cp != 0 && cp->configPipe(args)){
				//Run the pipe function (wrapper)
				*wD = cp->runPipeWrapper(*wD);
			} else {
				cout << "Failed to configure pipeline: " << args["pipeline"] << endl;
			}
		}
	}
	//If the pipeline was undefined...
	else {
		cout << "Failed to find a suitable pipeline, exiting..." << endl;
		return;
	}
	
	//Output the data using writeOutput library
	pipe = args.find("outputFile");
	if(pipe != args.end()){
		if (args["outputFile"] == "console"){
			//ws->writeConsole(wD);
		}
		
		ws->writeStats(wD->stats, args["outputFile"]);
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
			cout << "Failed to configure pipeline: " << args["pipeline"] << endl;
		}
	}
	
	args["pipeline"] = "distMatrix.neighGraph.rips.optPersistence";
	
	runPipeline(args, wD);
		
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
			cout << "Failed to configure pipeline: " << args["pipeline"] << endl;
		}
	}
	
	args["pipeline"] = "distMatrix.neighGraph.rips.optPersistence";
	
	do{
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
		
	return;
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
    
	//Create a pipePacket (datatype) to store the complex and pass between engines
    auto *wD = new pipePacket(args["complexType"], stod(args["epsilon"]), stoi(args["dimensions"]));	//wD (workingData)
	
	//Read data from inputFile CSV
    wD->workData.originalData = rs->readCSV(args["inputFile"]);

	//If data was found in the inputFile
	if(wD->workData.originalData.size() > 0){
		
		//Add data to our pipePacket
		wD->workData.originalData = wD->workData.originalData;
		
		auto ar = args["upscale"];
		if(ar == "true"){
			processUpscaleWrapper(args, wD);
		} else {
			processDataWrapper(args, wD);
		}
	} else {
		ap->printUsage();
	}
	
	
	
    return 0;
}
