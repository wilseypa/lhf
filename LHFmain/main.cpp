#include "mpi.h"
#include "omp.h"
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
#include <string> 

int nprocs,id;
 
struct sortBettis{
	bool operator()(bettiBoundaryTableEntry lhs, bettiBoundaryTableEntry rhs){
		return lhs.bettiDim < rhs.bettiDim || (lhs.bettiDim == rhs.bettiDim && lhs.death < rhs.death);
	}
};

void outputBettis(std::map<std::string, std::string> args, pipePacket* wD){
	auto *ws = new writeOutput();

	//Output the data using writeOutput library
	auto pipe = args.find("outputFile");
	if(pipe != args.end()){
		if (args["outputFile"] == "console"){
			//ws->writeConsole(wD);
		} else {
			ws->writeStats(wD->stats, args["outputFile"]);
			ws->writeBarcodes(wD->bettiTable, args["outputFile"]);	
		}
	}
	
	delete ws;
}

void runPipeline(std::map<std::string, std::string> args, pipePacket* wD){
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
			auto *cp = basePipe::newPipe(curFunct, args["complexType"]);
		
			//Check if the pipe was created and configure
			if(cp != 0 && cp->configPipe(args)){
				//Run the pipe function (wrapper)
				*wD = cp->runPipeWrapper(*wD);
			} else {
				std::cout<< cp <<std::endl;
				std::cout << "LHF runPipeline: Failed to configure pipeline: " << args["pipeline"] << std::endl;
			}
			
			delete cp;
		}
	}
	//If the pipeline was undefined...
	else {
		std::cout << "LHF runPipeline: Failed to find a suitable pipeline, exiting..." << std::endl;
		return;
	}
	
	//Output the data using writeOutput library
	outputBettis(args, wD);
}



void processDataWrapper(std::map<std::string, std::string> args, pipePacket* wD){
	
	//Start with the preprocessing function, if enabled
	auto pre = args["preprocessor"];
	if(pre != ""){
		auto *prePipe = preprocessor::newPreprocessor(pre);
		
		if(prePipe != 0 && prePipe->configPreprocessor(args)){
			*wD = prePipe->runPreprocessorWrapper(*wD);
		} else {
			std::cout << "LHF processData: Failed to configure pipeline: " << args["pipeline"] << std::endl;
		}
		
		delete prePipe;
	}
}	

std::vector<bettiBoundaryTableEntry> processIterUpscale(std::map<std::string, std::string> args, pipePacket* wD){
	
	//This function is called when the number of points in a partition are greater than the point threshold
	//	If the number of points in the new partitions are under the point threshold, continue
	//	If the number of points in any new partition is above the point threshold, recurse
	//
	//We may also want to limit the number of recursions for bigger data sets
	
	
	//Break this down into better code sections:
	
	
	//1. Local objects and storage
	
	//		Parameters	
	auto threshold = std::atoi(args["threshold"].c_str());
	auto maxEpsilon = std::atof(args["epsilon"].c_str());
	auto scalar = std::atof(args["scalar"].c_str());
	auto threads = std::atoi(args["threads"].c_str());
	auto clusters = std::atoi(args["clusters"].c_str());
	
	//		Referenced Libraries
	utils ut;
	
	//		Local Storage
	std::vector<bettiBoundaryTableEntry> mergedBettiTable;
	std::vector<bettiBoundaryTableEntry> partBettiTable[threads];
	auto originalDataSize = wD->originalData.size();
	
	//		Initalize a copy of the pipePacket
	auto iterwD = new pipePacket(args, args["complexType"]);
	iterwD->originalData = wD->originalData;
	iterwD->fullData = wD->fullData;
	
	//2. Partition the source point cloud separate datasets accordingly
	
	//		Run preprocessor pipeline to partition
	processDataWrapper(args, iterwD);
	
	//		Compute partition statistics for fuzzy partition distance
	auto maxRadius = ut.computeMaxRadius(clusters, iterwD->originalData, iterwD->fullData, iterwD->originalLabels);
	auto avgRadius = ut.computeAvgRadius(clusters, iterwD->originalData, iterwD->fullData, iterwD->originalLabels);
	
	//		Count the size of each partition for identifying source partitions when looking at the betti table results
	std::vector<unsigned> binCounts;
	for(unsigned a = 0; a < clusters; a++){
		binCounts.push_back(std::count(iterwD->originalLabels.begin(), iterwD->originalLabels.end(), a));
	}
	std::cout << "Bin Counts: ";
	ut.print1DVector(binCounts);
	
	//		Sort our fuzzy partitions into individual vectors
	args["scalarV"] = std::to_string(scalar*maxRadius);
	auto partitionedData = ut.separatePartitions(std::atof(args["scalarV"].c_str()), iterwD->originalData, iterwD->fullData, iterwD->originalLabels);
	std::cout << "Using scalar value: " << args["scalarV"] << std::endl;
	std::cout << "Partitions: " << partitionedData.second.size() << std::endl << "Partition Bin Counts: ";
	
	//		Append the centroid dataset to run in parallel as well
	partitionedData.second.push_back(iterwD->originalData);
	

	//		Get sizes of the new fuzzy partitions
	std::vector<unsigned> partitionsize;
	for(auto a: partitionedData.second)
		partitionsize.push_back(a.size());
	ut.print1DVector(partitionsize);
	
	//3. Process each partition using OpenMP to handle multithreaded scheduling
	
	std::cout << "Running with " << threads << " threads" << std::endl;
	
	//		This begins the parallel region; # threads are spun up within these braces
	#pragma omp parallel num_threads(threads)
	{
		
		//		Get the thread number with OpenMP - prevents simultaneous accesses to the betti tables
		int np = omp_get_thread_num();
		
		//		Schedule each thread to compute one of the partitions in any order
		//			When finished, thread will grab next iteration available
		#pragma omp for schedule(dynamic)
		for(int z = partitionedData.second.size()-1; z >= 0; z--){
			
			//Check if we are running the centroid dataset (there's no label associated)
			if(z == partitionedData.first.size()){
				//		Run on full dataset
				//		Set the pipeline to include the boundary and upscaling steps
				auto centArgs = args;
				if(args["upscale"] == "true" || args["upscale"] == "1")
					centArgs["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence.boundary.upscale";
				else
					centArgs["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence";
					
				if(partitionedData.second[z].size() > 0){
					wD->originalData = partitionedData.second[z];
					runPipeline(centArgs, wD);
					
					wD->complex->clear();
				} else 
					std::cout << "skipping full data, no centroids" << std::endl;
					
				//Determine if we need to upscale any additional boundaries based on the output of the centroid approximated PH
					
			} else if(partitionedData.second[z].size() > 0){
				
				//		Clone the pipePacket to prevent shared memory race conditions
				auto curwD = new pipePacket(args, args["complexType"]);	
				curwD->originalData = partitionedData.second[z];
				curwD->fullData = partitionedData.second[z];
				
				//		If the current partition is smaller than the threshold, process
				//			Otherwise recurse to reduce the number of points
				if(args["mode"] == "reduced" || partitionedData.second[z].size() < threshold){
					runPipeline(args, curwD);
				} else {
					curwD->bettiTable = processIterUpscale(args, curwD);
				}
				
				//4. Process and merge bettis - whether they are from runPipeline or IterUpscale
				bool foundExt = false;
				std::vector<bettiBoundaryTableEntry> temp;
				//ut.extractBoundaryPoints(curwD->bettiTable);

				//Remap the boundary indices into the original point space
				curwD->bettiTable = ut.mapPartitionIndexing(partitionedData.first[z], curwD->bettiTable);
              		
				for(auto betEntry : curwD->bettiTable){
					auto boundIter = betEntry.boundaryPoints.begin();
					
					//REWRITE::
					
					//The new (improved) approach to merging d0 bettis (pretty sure this works....)
					//	1. Use a binary array to track each point within the original partition
					//	2. Iterate the betti entries by weight, increasing
					//		a. If both indices are less than the partition size, check the binary array
					//			-If binary array for either of the two indices isn't filled, insert and fill all
					//		b. If one index is less than the partition size, the other greater, and this is the first instance of this
					//			-Add this to the connection list; this is the minimum connection outside of the partition
					//		c. If neither of the indices are less than the partition size, remove
					//	3. Once all entries have been iterated - if (b) was traversed there is a connection outside to another partition
					//		-If (b) was not traversed, need to add a {0, maxEps} entry for the independent component (Check this?)
					if(betEntry.boundaryPoints.size() > 0 && iterwD->originalLabels[*boundIter] == z){
						if(betEntry.bettiDim == 0){
							boundIter++;
							
							//Check if second entry is in the partition
							if(iterwD->originalLabels[*boundIter] == z){
								temp.push_back(betEntry);
							} else if(!foundExt){
								foundExt = true;
								temp.push_back(betEntry);
							}						
						} else{
							temp.push_back(betEntry);
						}
					}
				}
            	// //If we never found an external connection, add the infinite connection here
				// if(!foundExt){
				// 	bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {}, {} };
				// 	temp.push_back(des);
				// }
		
				for(auto newEntry : temp){
					bool found = false;
					for(auto curEntry : partBettiTable[np]){
						if(newEntry.death == curEntry.death && newEntry.boundaryPoints == curEntry.boundaryPoints){
							found = true;
						}
					}
					if(!found)
						partBettiTable[np].push_back(newEntry);
				}
				
				curwD->complex->clear();
				delete curwD->complex;
				delete curwD;
				
			} else 
				std::cout << "skipping" << std::endl;
		}
	}
	
	//5. Merge partitions, compute PH on original centroid dataset, and report results
	
	//		Merge partitioned betti tables together
	for(auto partTable : partBettiTable)
		mergedBettiTable.insert(mergedBettiTable.end(), partTable.begin(), partTable.end());	
		
	//		Merge bettis from the centroid based data
	for(auto betEntry : wD->bettiTable){
		if(betEntry.bettiDim > 0 ){
			mergedBettiTable.push_back(betEntry);
		}
	}
		
	// //		Add open d0 intervals for the remaining d0 bettis
	// auto addlIntervals = std::count_if(mergedBettiTable.begin(), mergedBettiTable.end(), [&](bettiBoundaryTableEntry const &i) { return ( i.bettiDim == 0); });
	// for(auto i = 0; i < originalDataSize - addlIntervals; i++)
	bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {}, {} };
	mergedBettiTable.push_back(des);

	iterwD->complex->clear();
	delete iterwD->complex;
	delete iterwD;
   // Return the final merged betti table for this iteration
	return mergedBettiTable;
}


void processUpscaleWrapper(std::map<std::string, std::string> args, pipePacket* wD){
	auto maxEpsilon = std::atof(args["epsilon"].c_str());
	auto scalar = std::atof(args["scalar"].c_str());
	//Start with the preprocessing function, if enabled
	auto *ws = new writeOutput();
	std::vector<bettiBoundaryTableEntry> mergedBettiTable;
	std::vector<bettiBoundaryTableEntry> finalMergedBettiTable;
	utils ut;
    
	// check if 1 process only
	if(nprocs==1){
		//Read our input data
		auto *rs = new readInput();
		wD->originalData = rs->readCSV(args["inputFile"]);
		wD->fullData = wD->originalData;

		auto *ws = new writeOutput();
		utils ut;
		
		args["mode"] = "reduced";
		auto mergedBettiTable = processIterUpscale(args, wD);
		// processReducedWrapper(args,wD);

		std::cout << std::endl << "_______Merged BETTIS_______" << std::endl;
		
		for(auto a : mergedBettiTable){
			std::cout << a.bettiDim << ",\t" << a.birth << ",\t" << a.death << ",\t";
			ut.print1DVector(a.boundaryPoints);
		}

		//Output the data using writeOutput library
		
		auto pipe = args.find("outputFile");
		if(pipe != args.end()){
			if (args["outputFile"] == "console"){
				//ws->writeConsole(wD);
			} else {
				ws->writeStats(wD->stats, args["outputFile"]);
				ws->writeBarcodes(mergedBettiTable, args["outputFile"]);		
			}
		}
		
		delete ws;
		return;
	}
	
	unsigned dimension;	
	std::vector<std::vector<double>> centroids;	
	std::vector<unsigned> binCounts;
	int minPartitions=0;
	int firstk=0;
	std::vector<unsigned> partitionsize;
	unsigned partitionsize_size=0;
	unsigned binCounts_size=0;
	std::vector<double> sdata;
	
	std::vector<unsigned> originalLabels;
	unsigned originalLabels_size;
	
	std::vector<unsigned> sdatalabel;
	int scountslabel[nprocs];
	int displslabel[nprocs];
	int bettiTableSize=0;
	unsigned boundary_size =0;
	int scounts[nprocs];
	int displs[nprocs];
	int maxsize=0;
	int maxsizelabel=0;
	double *receivedData;
	unsigned *receivedDataLabel;

	std::vector<unsigned> betti_dim;
	std::vector<double> betti_birth;
	std::vector<double> betti_death;
	std::vector<unsigned> betti_boundarysize;
	std::vector<unsigned> betti_boundaries;

	std::vector<bettiBoundaryTableEntry> mergedCentroidsBettiTable;           

	//Check if we are the master process for upscaling   
	if(id == 0){
			//Read our input data
		auto *rs = new readInput();
		wD->originalData = rs->readCSV(args["inputFile"]);
		wD->fullData = wD->originalData;
     
		//Partition the data with the configured preprocessor

		processDataWrapper(args, wD);
		
		originalLabels = wD->originalLabels;
		originalLabels_size = wD->originalLabels.size();
		
		//Separate our partitions for distribution
		auto maxRadius = ut.computeMaxRadius(std::atoi(args["clusters"].c_str()), wD->originalData, wD->fullData, wD->originalLabels);
		auto avgRadius = ut.computeAvgRadius(std::atoi(args["clusters"].c_str()), wD->originalData, wD->fullData, wD->originalLabels);
		centroids = wD->originalData;
		
             
		//Store the count of original points in each partition for merging
//		std::transform(std::begin(wD->originalLabels),std::end(wD->originalLabels),std::begin(wD->originalLabels),[maxlabel](unsigned x){return x+maxlabel;});
		for(unsigned a = 0; a < std::atoi(args["clusters"].c_str()); a++){
			binCounts.push_back(std::count(wD->originalLabels.begin(), wD->originalLabels.end(), a));
		}
       	auto partitionedData1 = ut.separatePartitions(scalar*maxRadius, wD->originalData, wD->fullData, wD->originalLabels);
        
        partitionedData1.second.push_back(centroids);
		//	Each node/slave will process at least 1 partition
		//	NOTE: the partition may contain points outside partition that are within 2*Rmax
		minPartitions = partitionedData1.second.size() / nprocs;
		firstk = partitionedData1.second.size() - (minPartitions*nprocs);
				
		dimension = partitionedData1.second[0][0].size();
		//Get the partition sizes
		for(auto a: partitionedData1.second)
			 partitionsize.push_back(a.size());
		//serializing all the data 
		for(auto a : partitionedData1.second)
			for(auto b : a)
				for(auto c : b)
					sdata.push_back(c);
	    //serialize all the labels
		for(auto a : partitionedData1.first)
			for(auto b : a)
				sdatalabel.push_back(b);
	    //intitalizing to 0
		for(int i=0;i<nprocs;i++){
			scounts[i] = 0;
			displs[i]= 0;
			scountslabel[i] = 0;
			displslabel[i]= 0;
		}
		
		//distributing partitions in scounts and scountslabel
		int temp=0,k=0;
         	for(auto a : partitionsize){
			if(k<firstk){
				if(temp>=minPartitions+1){
					temp=0;
					k++;
				}
			}
			else{
				if(temp>=minPartitions){
					temp=0;
					k++;
					}
			}
			scounts[k] += a*dimension;
			scountslabel[k] +=a;
			temp++;
		}
		k=0;
		maxsize=0;
		maxsizelabel=0;
		// find the maximum buffer size required for partition
        	for(auto  a : scounts){
			displs[k] = a + displs[k];
			k++;
			if(a>maxsize)
				maxsize = a;
		}
		k=0;
		// find maximum buffer size required by label partition	
		for(auto  a : scountslabel){
			displslabel[k] = a + displslabel[k];
			k++;
			if(a>maxsizelabel)
				maxsizelabel = a;
		}		
		partitionsize_size = partitionsize.size();
		binCounts_size = binCounts.size();
		std::cout << "partitionsize:";
		ut.print1DVector(partitionsize);
	}

	//	broadcasting all the required information by slaves.

	MPI_Bcast(&dimension,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);				//Data dimension
	MPI_Bcast(&partitionsize_size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);		//Partition size array size
	partitionsize.resize(partitionsize_size);
	MPI_Bcast(&partitionsize[0],partitionsize_size,MPI_UNSIGNED,0,MPI_COMM_WORLD); //Partition size array
	MPI_Bcast(&binCounts_size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);		//binCounts array size
	binCounts.resize(binCounts_size);
	MPI_Bcast(&binCounts[0],binCounts_size,MPI_UNSIGNED,0,MPI_COMM_WORLD); //binCounts size array
	
	MPI_Bcast(&originalLabels_size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);		//binCounts array size
	originalLabels.resize(originalLabels_size);
	MPI_Bcast(&originalLabels[0],originalLabels_size,MPI_UNSIGNED,0,MPI_COMM_WORLD); //binCounts size array
	
	MPI_Bcast(&firstk,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(&minPartitions,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&maxsize,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&scounts,nprocs,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&scountslabel,nprocs,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&maxsizelabel,1,MPI_INT,0,MPI_COMM_WORLD);

	if(id < firstk)
		minPartitions = minPartitions + 1;
	
    receivedData = (double *)malloc(maxsize*sizeof(double));
	receivedDataLabel = (unsigned *)malloc(maxsizelabel*sizeof(unsigned));
	for(int i=0;i<maxsize;i++){
			receivedData[i] = 0.0;
	}
	for(int i=0;i<maxsizelabel;i++){
			receivedDataLabel[i] = 0;
	}

	//	scatter the partitions to slaves.
  MPI_Scatterv( &sdata[0],scounts,displs, MPI_DOUBLE, receivedData, maxsize, MPI_DOUBLE,0, MPI_COMM_WORLD);
  MPI_Scatterv( &sdatalabel[0],scountslabel,displslabel, MPI_UNSIGNED, receivedDataLabel, maxsizelabel, MPI_UNSIGNED,0, MPI_COMM_WORLD);
  
    if(minPartitions>0){
		//NOTE: need to have dynamic partition size; whether that means serializing and sending
		//	the partition table and size or dynamically allocating the partitionsize vector here (push_back)
	//		std::cout<<dimension<<" "<<minPartitions<<" "<<id<<std::endl;
		int i=0;
		int displacement=0;
		if(id < firstk)
			displacement = id*minPartitions;
		else
			displacement = firstk*(minPartitions+1) + (id-firstk)*minPartitions;
		
		std::vector<std::vector<std::vector<double>>> partitionData;
		std::vector<std::vector<unsigned>> labels;
		int index=0;
		int p=0;
		for(int i=0;i<minPartitions;i++){
			std::vector<std::vector<double>> partition;
			for(int j=0;j<partitionsize[displacement+i];j++){
				std::vector<double> row;
				for(int k=0;k<dimension;k++){
				   row.push_back(receivedData[p++]);
				}
			partition.push_back(row);	
			}
			partitionData.push_back(partition);
		//	std::cout<<"data = "<<partition.size()<<" "<<"id ="<<id<<std::endl;			
		}
        
		
		p=0;
                for(int i=0;i<minPartitions;i++){
			std::vector<unsigned> partitionlabel;
			for(int j=0;j<partitionsize[displacement+i];j++){
				partitionlabel.push_back(receivedDataLabel[p++]);
			}
			labels.push_back(partitionlabel);
		//	std::cout<<"label ="<<partitionlabel.size()<<" "<<"id ="<<id<<std::endl;
		}
	auto partitionedData = make_pair(labels,partitionData);	
			
	auto threshold = std::atoi(args["threshold"].c_str());
	auto maxEpsilon = std::atof(args["epsilon"].c_str());
	auto scalar = std::atof(args["scalar"].c_str());
	auto threads = std::atoi(args["threads"].c_str());
	auto clusters = std::atoi(args["clusters"].c_str());
	
	//		Referenced Libraries
	utils ut;
	//		Local Storage
		//		Initalize a copy of the pipePacket
	auto iterwD =new pipePacket(args, args["complexType"]);
	iterwD->originalData = wD->originalData;
	iterwD->fullData = wD->fullData;
	
	std::vector<bettiBoundaryTableEntry> mergedBettiTable;
	std::vector<bettiBoundaryTableEntry> partBettiTable[threads];
	auto originalDataSize = wD->originalData.size();
		
	//3. Process each partition using OpenMP to handle multithreaded scheduling
	
	std::cout << "Running with " << threads << " threads" << std::endl;
	
	//		This begins the parallel region; # threads are spun up within these braces
	#pragma omp parallel num_threads(threads)
	{
		
		//		Get the thread number with OpenMP - prevents simultaneous accesses to the betti tables
		int np = omp_get_thread_num();
		
		//		Schedule each thread to compute one of the partitions in any order
		//			When finished, thread will grab next iteration available
		#pragma omp for schedule(dynamic)
		for(int z = partitionedData.second.size()-1; z >= 0; z--){
			
			//Check if we are running the centroid dataset (there's no label associated)
			if(partitionedData.second[z].size() == clusters){
				
				std::cout<<"Running centroid "<<clusters<<"  id = "<<id<<std::endl;
				//		Run on full dataset
				//		Set the pipeline to include the boundary and upscaling steps
				auto centArgs = args;
				if(args["upscale"] == "true" || args["upscale"] == "1")
					centArgs["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence.boundary.upscale";
				else
					centArgs["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence";
					
				//Run against the original dataset
			
				if(partitionedData.second[z].size() > 0){
					iterwD->originalData = partitionedData.second[z];
					iterwD->originalData = partitionedData.second[z];
					runPipeline(centArgs, iterwD);
					delete iterwD->complex;
				} else 
					std::cout << "skipping full data, no centroids" << std::endl;
					
				//Determine if we need to upscale any additional boundaries based on the output of the centroid approximated PH
					
			} else if(partitionedData.second[z].size() > 0){
				
				//		Clone the pipePacket to prevent shared memory race conditions
			
				auto curwD = new  pipePacket(args, args["complexType"]);	
				curwD->originalData = partitionedData.second[z];
				curwD->fullData = partitionedData.second[z];
				
				//		If the current partition is smaller than the threshold, process
				//			Otherwise recurse to reduce the number of points
				if(args["mode"] == "reduced" || partitionedData.second[z].size() < threshold){
					runPipeline(args, curwD);
				} else {
					curwD->bettiTable = processIterUpscale(args, curwD);
     			}
				delete curwD->complex;
				//4. Process and merge bettis - whether they are from runPipeline or IterUpscale
				bool foundExt = false;
				std::vector<bettiBoundaryTableEntry> temp;
				//ut.extractBoundaryPoints(curwD.bettiTable);

				//Remap the boundary indices into the original point space
				curwD->bettiTable = ut.mapPartitionIndexing(partitionedData.first[z], curwD->bettiTable);
               
				//ut.print1DVector(originalLabels);
				
				for(auto betEntry : curwD->bettiTable){
					auto boundIter = betEntry.boundaryPoints.begin();
					
					//REWRITE::
					
					//The new (improved) approach to merging d0 bettis (pretty sure this works....)
					//	1. Use a binary array to track each point within the original partition
					//	2. Iterate the betti entries by weight, increasing
					//		a. If both indices are less than the partition size, check the binary array
					//			-If binary array for either of the two indices isn't filled, insert and fill all
					//		b. If one index is less than the partition size, the other greater, and this is the first instance of this
					//			-Add this to the connection list; this is the minimum connection outside of the partition
					//		c. If neither of the indices are less than the partition size, remove
					//	3. Once all entries have been iterated - if (b) was traversed there is a connection outside to another partition
					//		-If (b) was not traversed, need to add a {0, maxEps} entry for the independent component (Check this?)
					if(betEntry.boundaryPoints.size() > 0 && originalLabels[*boundIter] == displacement+z){
						if(betEntry.bettiDim == 0){
							boundIter++;
							
							//Check if second entry is in the partition
							if(originalLabels[*boundIter] == displacement+z){
								temp.push_back(betEntry);
							} else if(!foundExt){
								foundExt = true;
								temp.push_back(betEntry);
							}						
						} else{
							temp.push_back(betEntry);
						}
					}
				}

				// //If we never found an external connection, add the infinite connection here
				// if(!foundExt){
				// 	bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {}, {} };
				// 	temp.push_back(des);
				// }
		      	for(auto newEntry : temp){
					bool found = false;
					for(auto curEntry : partBettiTable[np]){
						if(newEntry.death == curEntry.death && newEntry.boundaryPoints == curEntry.boundaryPoints){
							found = true;
						}
					}
					if(!found)
						partBettiTable[np].push_back(newEntry);
				}
								
			} else 
				std::cout << "skipping" << std::endl;
		}
	}
		
	//		Merge partitioned betti tables together
	for(auto partTable : partBettiTable)
		mergedBettiTable.insert(mergedBettiTable.end(), partTable.begin(), partTable.end());	
		
		
	//		Merge bettis from the centroid based data
	if(id==nprocs-1)
	{
	for(auto betEntry : iterwD->bettiTable){
	//	if(betEntry.bettiDim > 0 ){
			mergedBettiTable.push_back(betEntry);
	//	}
	}
	}
	
	// //		Add open d0 intervals for the remaining d0 bettis
	// auto addlIntervals = std::count_if(mergedBettiTable.begin(), mergedBettiTable.end(), [&](bettiBoundaryTableEntry const &i) { return ( i.bettiDim == 0); });
	// for(auto i = 0; i < originalDataSize - addlIntervals; i++)
	bettiBoundaryTableEntry des = { 0, 0, maxEpsilon, {} };
	mergedBettiTable.push_back(des);
		
		for(auto bet : mergedBettiTable){
			betti_dim.push_back(bet.bettiDim);
			betti_birth.push_back(bet.birth);
			betti_death.push_back(bet.death);
			betti_boundarysize.push_back(bet.boundaryPoints.size());
			betti_boundaries.insert(betti_boundaries.end(),bet.boundaryPoints.begin(),bet.boundaryPoints.end());
	  
		}
		
		bettiTableSize = betti_dim.size();
		boundary_size = 0;
		for(auto a : betti_boundarysize)
			boundary_size +=a;
	}	
    else{
		std::cout<<"**********************\n"<<std::endl;
	}
	//gether all the betties from slaves
	int totalsize =0;
	MPI_Reduce(&bettiTableSize, &totalsize, 1, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);
	std::vector<int> num(nprocs);
	MPI_Gather( &bettiTableSize, 1, MPI_INT, &num[0], 1, MPI_INT,0, MPI_COMM_WORLD); 
    std::vector<int> displsg(nprocs);
	displsg[0]=0;

	for(int i=1;i<=nprocs;i++)
		displsg[i] = displsg[i-1] + num[i-1];

	std::vector<unsigned> recvbuffdim(totalsize);
	std::vector<double> recvbuffbirth(totalsize);
	std::vector<double> recvbuffdeath(totalsize);
	std::vector<unsigned> recvbuffboundaries_size(totalsize);
	MPI_Gatherv(&betti_dim[0],bettiTableSize, MPI_UNSIGNED, &recvbuffdim[0], &num[0], &displsg[0],MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Gatherv(&betti_birth[0],bettiTableSize, MPI_DOUBLE, &recvbuffbirth[0], &num[0], &displsg[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gatherv(&betti_death[0],bettiTableSize, MPI_DOUBLE, &recvbuffdeath[0], &num[0], &displsg[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gatherv(&betti_boundarysize[0],bettiTableSize, MPI_UNSIGNED, &recvbuffboundaries_size[0], &num[0], &displsg[0],MPI_UNSIGNED,0,MPI_COMM_WORLD);
    int totalboundarysize =0;
	MPI_Reduce(&boundary_size, &totalboundarysize, 1, MPI_UNSIGNED, MPI_SUM, 0,MPI_COMM_WORLD);
	
	std::vector<unsigned> recvbuffboundaries(totalboundarysize);
	std::vector<int> num1(nprocs);
	MPI_Gather(&boundary_size, 1, MPI_UNSIGNED, &num1[0], 1, MPI_UNSIGNED,0, MPI_COMM_WORLD);
  
    std::vector<int> displsgb(nprocs);
	displsgb[0]=0;

	for(int i=1;i<=nprocs;i++){
		displsgb[i] = displsgb[i-1] + num1[i-1];
	}
	MPI_Gatherv(&betti_boundaries[0],boundary_size, MPI_UNSIGNED, &recvbuffboundaries[0], &num1[0], &displsgb[0],MPI_UNSIGNED,0,MPI_COMM_WORLD);

	if(id==0){
		//master prune out the duplicates across slaves.
		
		int beg=0;
		for(int i=0;i<totalsize;i++){
			bettiBoundaryTableEntry bettiEntry;
			bettiEntry.bettiDim = recvbuffdim[i];
			bettiEntry.birth = recvbuffbirth[i];
			bettiEntry.death = recvbuffdeath[i];
			for(int bi = beg;bi<(beg+recvbuffboundaries_size[i]);bi++)
				bettiEntry.boundaryPoints.insert(recvbuffboundaries[bi]);
			
			beg +=recvbuffboundaries_size[i];
			bool found = false;
			for(auto curEntry : finalMergedBettiTable){
				if(bettiEntry.death == curEntry.death && bettiEntry.boundaryPoints == curEntry.boundaryPoints){
					found = true;
					}
			}
			if(!found)
				finalMergedBettiTable.push_back(bettiEntry);
		}


		sort(finalMergedBettiTable.begin(), finalMergedBettiTable.end(), sortBettis());	
		std::cout << std::endl << "_______Merged BETTIS_______" << std::endl;
		
		for(auto a : finalMergedBettiTable){
			std::cout << a.bettiDim << ",\t" << a.birth << ",\t" << a.death << ",\t";
			ut.print1DVector(a.boundaryPoints);
		}
			
		//Output the data using writeOutput library
		auto pipe = args.find("outputFile");
		if(pipe != args.end()){
			if (args["outputFile"] == "console"){
				//ws->writeConsole(wD);
			} else {
				ws->writeStats(wD->stats, args["outputFile"]);
				ws->writeBarcodes(finalMergedBettiTable, args["outputFile"]);
				
			}
		}		
		
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
	
	auto *rs = new readInput();
	auto *ap = new argParser();
	
	//Parse the command-line arguments
	auto args = ap->parse(argc, argv);
	
	//Determine what pipe we will be running
	ap->setPipeline(args);
	
	//Create a pipePacket (datatype) to store the complex and pass between engines
	auto *wD = new pipePacket(args, args["complexType"]);	//wD (workingData)
	
	if(args["pipeline"] != "slidingwindow" && args["pipeline"] != "naivewindow" && args["mode"] != "mpi"){
		//Read data from inputFile CSV
		wD->originalData = rs->readCSV(args["inputFile"]);
		wD->fullData = wD->originalData;
	}
	
	//If data was found in the inputFile
	if(wD->originalData.size() > 0 || args["pipeline"] == "slidingwindow" || args["pipeline"] == "naivewindow" || args["mode"] == "mpi"){

		//Add data to our pipePacket
		wD->originalData = wD->originalData;
		
		if(args["mode"] == "mpi"){
			
			MPI_Init(&argc,&argv);
			MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
			MPI_Comm_rank(MPI_COMM_WORLD,&id);
			
			processUpscaleWrapper(args, wD);
			
			MPI_Finalize();

		} else if(args["mode"] == "reduced" || args["mode"] == "iterUpscale" || args["mode"] == "iter"){	
			wD->bettiTable = processIterUpscale(args,wD);
			sort(wD->bettiTable.begin(), wD->bettiTable.end(), sortBettis());
			utils ut;
			
			if(args["debug"] == "1" || args["debug"] == "true"){
				std::cout << std::endl << "_______Merged BETTIS_______" << std::endl;
		
				for(auto a : wD->bettiTable){
					std::cout << a.bettiDim << ",\t" << a.birth << ",\t" << a.death << ",\t";
					ut.print1DVector(a.boundaryPoints);
				}
			}
			
			//Output the data using writeOutput library
			outputBettis(args, wD);
			
		} else {
			processDataWrapper(args, wD);
			runPipeline(args, wD);
		}
	} else {
		ap->printUsage();
	}
	
	wD->complex->clear();
	delete wD->complex;
	delete rs, delete ap, delete wD;

	return 0;
}
