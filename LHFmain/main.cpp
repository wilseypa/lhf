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

int nprocs,id;
 
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
	std::vector<bettiBoundaryTableEntry> mergedBettiTable;
    utils ut;
    
    // check if 1 process only
      if(nprocs==1)
      {
		 //Read our input data
		auto *rs = new readInput();
		wD->originalData = rs->readCSV(args["inputFile"]);
		wD->fullData = wD->originalData;
		
		processReducedWrapper(args,wD);
		 
		return ;
		
		  }
    //Check if we are the master process for upscaling   
 	if(id == 0){
		std::vector<bettiBoundaryTableEntry> mergedMasterBettiTable;
		
		//Read our input data
		auto *rs = new readInput();
		wD->originalData = rs->readCSV(args["inputFile"]);
		wD->fullData = wD->originalData;

		//Partition the data with the configured preprocessor
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
		//Separate our partitions for distribution
		auto maxRadius = ut.computeMaxRadius(std::atoi(args["clusters"].c_str()), wD->originalData, wD->fullData, wD->originalLabels);
		auto centroids = wD->originalData;
		std::cout << "Using maxRadius: " << maxRadius << std::endl;
		
		auto partitionedData = ut.separatePartitions(2*maxRadius, wD->originalData, wD->fullData, wD->originalLabels);
		
		//	Each node/slave will process at least 1 partition
		//		NOTE: the partition may contain points outside partition that are within 2*Rmax
		int minPartitions = partitionedData.second.size() / (nprocs-1);
		int firstk = partitionedData.second.size() - (minPartitions*(nprocs-1));
				
		int dimension = partitionedData.second[0][0].size();
		std::vector<unsigned> partitionsize;
		
		//Get the partition sizes
		for(auto a: partitionedData.second)
			 partitionsize.push_back(a.size());

		int partitionsize_size = partitionsize.size();
		ut.print1DVector(partitionsize);
		MPI_Request req;
		//Sending the data dimensions to slaves	 
		for(int i=1;i<nprocs;i++){			
			MPI_Isend(&dimension,1,MPI_INT,i,1,MPI_COMM_WORLD,&req);
			MPI_Isend(&partitionsize_size,1,MPI_INT,i,1,MPI_COMM_WORLD,&req);
			MPI_Isend(&partitionsize[0],partitionsize_size,MPI_UNSIGNED,i,1,MPI_COMM_WORLD,&req);
			if(i<=firstk){
				int part_size = minPartitions +1;
				MPI_Isend(&part_size,1,MPI_INT,i,1,MPI_COMM_WORLD,&req);
			}
			else
				MPI_Isend(&minPartitions,1,MPI_INT,i,1,MPI_COMM_WORLD,&req);			
		}
		//Sending Partitions to slaves
		for(int k =1; k < nprocs; k++){
			for(int i=k-1; i<partitionsize_size; i += (nprocs-1)){
				std::vector<double> part = ut.serialize(partitionedData.second[i]);
				MPI_Isend(&part[0], part.size(), MPI_DOUBLE, k, 1, MPI_COMM_WORLD,&req);
				MPI_Isend(&partitionedData.first[i],partitionedData.first[i].size(),MPI_UNSIGNED,k,1,MPI_COMM_WORLD,&req);
			}
		}
		
		std::cout << "Full Data: " << centroids.size() << std::endl;
		if(centroids.size() > 0){
			std::cout << "Running Pipeline with : " << centroids.size() << " vectors" << std::endl;
			wD->originalData = centroids;
			runPipeline(args, wD);
		
			wD->complex->clear();
		} else 
			std::cout << "skipping" << std::endl;
		std::vector<unsigned> betti_dim;
		std::vector<double> betti_birth;
		std::vector<double> betti_death;
		std::vector<unsigned> betti_boundarysize;
		std::vector<unsigned> betti_boundaries;
		for(int k=1; k < nprocs; k++){
	 		int p=0;
			int sizeOfBettiTable;
			
			MPI_Status statusr1;
			MPI_Request reqr1;
			MPI_Irecv(&sizeOfBettiTable,1,MPI_INT,k,p++,MPI_COMM_WORLD,&reqr1);
			MPI_Wait(&reqr1,&statusr1);
			std::cout << "Receiving data from slaves..."<<k<< std::endl;
			
			betti_dim.resize(sizeOfBettiTable);
			betti_birth.resize(sizeOfBettiTable);
			betti_death.resize(sizeOfBettiTable);
			betti_boundarysize.resize(sizeOfBettiTable);
		    
			MPI_Status statusr2;
			MPI_Request reqr2;
			MPI_Irecv(&betti_dim[0],sizeOfBettiTable,MPI_UNSIGNED,k,p++,MPI_COMM_WORLD,&reqr2);
			
			MPI_Status statusr3;
			MPI_Request reqr3;
			MPI_Irecv(&betti_birth[0],sizeOfBettiTable,MPI_DOUBLE,k,p++,MPI_COMM_WORLD,&reqr3);
			
			MPI_Status statusr4;
			MPI_Request reqr4;
    		MPI_Irecv(&betti_death[0],sizeOfBettiTable,MPI_DOUBLE,k,p++,MPI_COMM_WORLD,&reqr4);
			
			MPI_Status statusr5;
			MPI_Request reqr5;
			MPI_Irecv(&betti_boundarysize[0],sizeOfBettiTable,MPI_INT,k,p++,MPI_COMM_WORLD,&reqr5);
			
			MPI_Request reqr6;
			MPI_Status statusr6;
			int bound_size =0;
			MPI_Irecv(&bound_size,1,MPI_INT,k,p++,MPI_COMM_WORLD,&reqr6);
			MPI_Wait(&reqr6,&statusr6);
			
			betti_boundaries.resize(bound_size);
			MPI_Status statusr7;
			MPI_Request reqr7;
			MPI_Irecv(&betti_boundaries[0],bound_size,MPI_UNSIGNED,k,p++,MPI_COMM_WORLD,&reqr7);
			
			MPI_Wait(&reqr2,&statusr2);
			MPI_Wait(&reqr3,&statusr3);
			MPI_Wait(&reqr4,&statusr4);
			MPI_Wait(&reqr5,&statusr5);
			MPI_Wait(&reqr7,&statusr7);
			
			
			int beg = 0;
			for(int i=0;i<sizeOfBettiTable;i++){
				bettiBoundaryTableEntry bettiEntry;
				bettiEntry.bettiDim = betti_dim[i];
				bettiEntry.birth = betti_birth[i];
				bettiEntry.death = betti_death[i];
				
				for(int bi = beg;bi<(beg+betti_boundarysize[i]);bi++)
					bettiEntry.boundaryPoints.insert(betti_boundaries[bi]);
				
				beg +=betti_boundarysize[i];
				mergedMasterBettiTable.push_back(bettiEntry);
			}

		betti_dim.clear();
		betti_birth.clear();
		betti_death.clear();
		betti_boundarysize.clear();
		betti_boundaries.clear();


			/*
			for(int i=0;i<sizeOfBettiTable;i++){
				int boundarysize;		
				MPI_Irecv(&boundarysize,1,MPI_INT,k,p++,MPI_COMM_WORLD,&req);
				bettiBoundaryTableEntry bettiEntry;
				MPI_Irecv(&bettiEntry.bettiDim,1,MPI_UNSIGNED,k,p++,MPI_COMM_WORLD,&req);
				MPI_Irecv(&bettiEntry.birth,1,MPI_DOUBLE,k,p++,MPI_COMM_WORLD,&req);
				MPI_Irecv(&bettiEntry.death,1,MPI_DOUBLE,k,p++,MPI_COMM_WORLD,&req);
				std:vector<unsigned> boundaryvector(boundarysize);
				boundaryvector.resize(boundarysize);
				MPI_Irecv(&boundaryvector[0],boundarysize,MPI_UNSIGNED,k,p++,MPI_COMM_WORLD,&req);
				for(auto a:boundaryvector)
					bettiEntry.boundaryPoints.insert(a);
				
				mergedMasterBettiTable.push_back(bettiEntry);
			}
			*/
		}
		
	
		std::cout << std::endl << "_______BETTIS_______" << std::endl;
	
		//for(auto a : mergedMasterBettiTable){
			//std::cout << a.bettiDim << ",\t" << a.birth << ",\t" << a.death << ",\t";
			//ut.print1DVector(a.boundaryPoints);
		//}
		std::cout << std::endl << "mergedMasterBettiTable size is :: " << mergedMasterBettiTable.size() << std::endl;
			
			

	//SLAVE PROCESS:

	} else {
		int dim;
		//NOTE: need to have dynamic partition size; whether that means serializing and sending
		//	the partition table and size or dynamically allocating the partsize vector here (push_back)
		std::vector<unsigned> partsize;
		int partsize_size;
		MPI_Status status;
		MPI_Request req;
		//MPI_RECV( &Data, Size, 
		//Get the dimension of data coming in
		MPI_Status status1;
		MPI_Request req1;
		MPI_Irecv(&dim,1,MPI_INT,0,1,MPI_COMM_WORLD,&req1);
		
		//Get the total number of partitions
		MPI_Status status2;
		MPI_Request req2;
		MPI_Irecv(&partsize_size,1,MPI_INT,0,1,MPI_COMM_WORLD,&req2);
		MPI_Wait(&req2,&status2);
		partsize.resize(partsize_size);
		
		//Get the partition sizes
		MPI_Status status3;
		MPI_Request req3;
		MPI_Irecv(&partsize[0],partsize_size,MPI_UNSIGNED,0,1,MPI_COMM_WORLD,&req3);
		
		int no_of_partition;
		MPI_Status status4;
		MPI_Request req4;
		MPI_Irecv(&no_of_partition,1,MPI_INT,0,1,MPI_COMM_WORLD,&req4);
		
		std::vector<double> flatPartitions;
		std::vector<unsigned> labels;

		int p=0;	
    	MPI_Wait(&req1,&status1);
    	MPI_Wait(&req3,&status3);
		MPI_Wait(&req4,&status4);
		//std::vector<bettiBoundaryTableEntry> mergedBettiTableSlave;
		std::vector<unsigned> betti_dim;
		std::vector<double> betti_birth;
		std::vector<double> betti_death;
		std::vector<unsigned> betti_boundarysize;
		std::vector<unsigned> betti_boundaries;
		for(unsigned z = 0; z < no_of_partition; z++){
			int rsize = partsize[z*(nprocs-1)+(id-1)]*dim;
			flatPartitions.resize(rsize);
			MPI_Request req5;
			MPI_Status status5;
			MPI_Irecv(&flatPartitions[0],rsize,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&req5);
			MPI_Wait(&req5,&status5);	
			auto partitionedData = ut.deserialize(flatPartitions,dim);
			int lsize = partsize[z*(nprocs-1)+(id-1)];
			labels.resize(lsize);
			MPI_Request req6;
			MPI_Status status6;			
			MPI_Irecv(&labels[0],lsize,MPI_UNSIGNED,0,1,MPI_COMM_WORLD,&req6);
			MPI_Wait(&req6,&status6);
			
		if(partitionedData.size() > 0){
			std::cout << "Running Pipeline with : " << partitionedData.size() << " vectors" << " id :: "<<id<<std::endl;
			wD->originalData = partitionedData;
			runPipeline(args, wD);
			wD->complex->clear();
			//Map partitions back to original point indexing
			ut.mapPartitionIndexing(labels, wD->bettiTable);
			for(auto bet : wD->bettiTable){
				//mergedBettiTableSlave.push_back(bet);
		        betti_dim.push_back(bet.bettiDim);
				betti_birth.push_back(bet.birth);
				betti_death.push_back(bet.death);
				betti_boundarysize.push_back(bet.boundaryPoints.size());
				betti_boundaries.insert(betti_boundaries.end(),bet.boundaryPoints.begin(),bet.boundaryPoints.end());
			}
		} else 
			std::cout << "skipping" << std::endl;
    	}
		
			//Sending of merged betti results back to master (but this should be done after all partitions are finished)
			int bettiTableSize = betti_dim.size();
			MPI_Status statuss1;
			MPI_Request reqs1;
			MPI_Isend(&bettiTableSize,1,MPI_INT,0,p++,MPI_COMM_WORLD,&reqs1);
			MPI_Wait(&reqs1,&statuss1);
			
			MPI_Status statuss2;
			MPI_Request reqs2;
			MPI_Isend(&betti_dim[0],bettiTableSize,MPI_UNSIGNED,0,p++,MPI_COMM_WORLD,&reqs2);
			
			MPI_Status statuss3;
			MPI_Request reqs3;
			MPI_Isend(&betti_birth[0],bettiTableSize,MPI_DOUBLE,0,p++,MPI_COMM_WORLD,&reqs3);
			
			MPI_Status statuss4;
			MPI_Request reqs4;
			MPI_Isend(&betti_death[0],bettiTableSize,MPI_DOUBLE,0,p++,MPI_COMM_WORLD,&reqs4);
			
			MPI_Status statuss5;
			MPI_Request reqs5;
			MPI_Isend(&betti_boundarysize[0],bettiTableSize,MPI_UNSIGNED,0,p++,MPI_COMM_WORLD,&reqs5);
			
			int bound_size = 0;
			for(auto e :betti_boundarysize)
				bound_size += e;
			
			MPI_Status statuss6;
			MPI_Request reqs6;
			MPI_Isend(&bound_size,1,MPI_INT,0,p++,MPI_COMM_WORLD,&reqs6);
			MPI_Wait(&reqs6,&statuss6);
			
			MPI_Status statuss7;
			MPI_Request reqs7;
			MPI_Isend(&betti_boundaries[0],bound_size,MPI_UNSIGNED,0,p++,MPI_COMM_WORLD,&reqs7);
			MPI_Wait(&reqs2,&statuss2);
			MPI_Wait(&reqs3,&statuss3);
			MPI_Wait(&reqs4,&statuss4);
			MPI_Wait(&reqs5,&statuss5);
			MPI_Wait(&reqs7,&statuss7);
			
			
			/*
			for(auto betEntry : mergedBettiTableSlave){  
				int boundary_size = betEntry.boundaryPoints.size();	
				MPI_Isend(&boundary_size,1,MPI_INT,0,p++,MPI_COMM_WORLD,&req);
				MPI_Isend(&betEntry.bettiDim,1,MPI_UNSIGNED,0,p++,MPI_COMM_WORLD,&req);
				MPI_Isend(&betEntry.birth,1,MPI_DOUBLE,0,p++,MPI_COMM_WORLD,&req);
				MPI_Isend(&betEntry.death,1,MPI_DOUBLE,0,p++,MPI_COMM_WORLD,&req);
				vector<unsigned> boundaryvector;
				boundaryvector.assign(betEntry.boundaryPoints.begin(), betEntry.boundaryPoints.end());
				MPI_Isend(&boundaryvector[0],boundary_size,MPI_UNSIGNED,0,p++,MPI_COMM_WORLD,&req);
			}	
            */	    	
	}
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
    MPI_Init(&argc,&argv);
   	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	
	auto *rs = new readInput();
	auto *ap = new argParser();
    
    //Parse the command-line arguments
    auto args = ap->parse(argc, argv);
    
    //Determine what pipe we will be running
    ap->setPipeline(args);
    
 //   for(auto z : args)
	//	std::cout << z.first << "\t" << z.second << std::endl;
    
	//Create a pipePacket (datatype) to store the complex and pass between engines
    auto *wD = new pipePacket(args, args["complexType"]);	//wD (workingData)
	
	if(args["pipeline"] != "slidingwindow"&&args["mode"] != "mpi"){
		//Read data from inputFile CSV
		wD->originalData = rs->readCSV(args["inputFile"]);
		wD->fullData = wD->originalData;
	}
	
	//If data was found in the inputFile
	if(wD->originalData.size() > 0 || args["pipeline"] == "slidingwindow" || args["mode"] == "mpi"){
		
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
	
	
	 MPI_Finalize();
    return 0;
}
