/*
@file LHF.hpp
@brief Definition of LHF class, which runs the pipeline and outputs betti numbers.

*/
#include "mpi.h"
#include "LHF.hpp"
#include "omp.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <thread>
#include <string>
#include <numeric>



template<typename nodeType>
void LHF<nodeType>::outputBettis(std::map<std::string, std::string> args, pipePacket<nodeType> &wD){
    /**
        @brief Output betti numbers to a file or the console using the writeOutput library.

        @tparam nodeType The type of node in the data set.

        @param args A map containing the arguments for the pipeline.

        @param wD A pipePacket containing the output of the pipeline.
    */
    
	//Output the data using writeOutput library
	auto pipe = args.find("outputFile");
	if (pipe != args.end()){
		if (args["outputFile"] == "console"){
			writeOutput::writeConsole(wD.bettiTable);

			//Check if debug mode for runLog, stats
			pipe = args.find("debug");
			if (pipe != args.end() && args["debug"] == "1"){
				writeOutput::writeRunLog(wD.runLog, args["outputFile"]);
				writeOutput::writeStats(wD.stats, args["outputFile"]);
			}

		} else {
			sort(wD.bettiTable.begin(), wD.bettiTable.end(), sortBettis());
			writeOutput::writeBarcodes(wD.bettiTable, args["outputFile"]);

			//Check if debug mode for runLog, stats
			pipe = args.find("debug");
			if (pipe != args.end() && args["debug"] == "1"){
				writeOutput::writeRunLog(wD.runLog, args["outputFile"]);
				writeOutput::writeStats(wD.stats, args["outputFile"]);
			}
		}
	}
}


template<typename nodeType>
void LHF<nodeType>::runPipeline(std::map<std::string, std::string> args, pipePacket<nodeType>&wD){
    /**
        @brief Runs the pipeline with the specified arguments and data packet.

        The pipeline consists of a sequence of components connected by pipes

        that process the data packet in a specific order.

        The pipeline function names are separated by dots and passed as an argument.

        @tparam nodeType The data type of the nodes in the pipeline.

        @param args A map of arguments to configure the pipeline components.

        @param wD The input and output data packet for the pipeline.
    */
    
	// Begin processing parts of the pipeline
	// DataInput -> A -> B -> ... -> DataOutput
	// Parsed by "." -> i.e. A.B.C.D
	auto pipe = args.find("pipeline");
	
	if (pipe != args.end()){

		auto pipeFuncts = std::string(args["pipeline"]);
		auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;

		//Start the timer for time passed during the pipeline
		auto startTime = std::chrono::high_resolution_clock::now();

		//For each '.' separated pipeline function (count of '.' + 1 -> lim)
		for (unsigned i = 0; i < lim; i++){

			auto curFunct = pipeFuncts.substr(0, pipeFuncts.find('.'));
			pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);
			//Build the pipe component, configure and run
			auto cp = basePipe<nodeType>::newPipe(curFunct, args["complexType"]);

			//Check if the pipe was created and configure
			if (cp != 0 && cp->configPipe(args)){
				//Run the pipe function (wrapper)
				auto startTimeEach = std::chrono::high_resolution_clock::now();

             	cp->runPipeWrapper(wD);
             	auto endTimeEach = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> elapsedEach = endTimeEach - startTimeEach;
                std::cout<<" Pipe "<<curFunct<<"  Time Elapsed ::"<<std::to_string(elapsedEach.count() / 1000.0)<<std::endl;
			}
			else{
				std::cout << "LHF runPipeline: Failed to configure pipeline: " << args["pipeline"] << std::endl;
			}

			delete cp;
		}

		//Stop the timer for time passed during the pipeline
		auto endTime = std::chrono::high_resolution_clock::now();

		//Calculate the duration (physical time) for the pipeline
		std::chrono::duration<double, std::milli> elapsed = endTime - startTime;

		//Log the current execution to runLog
		wD.runLog += writeOutput::logRun(args, wD.ident, wD.getStats(), std::to_string(elapsed.count() / 1000.0));
	}
	//If the pipeline was undefined...
	else{
		std::cout << "LHF runPipeline: Failed to find a suitable pipeline, exiting..." << std::endl;
		return;
	}

	outputBettis(args, wD);
}

template<typename nodeType>
void LHF<nodeType>::runPreprocessor(std::map<std::string, std::string>& args, pipePacket<nodeType>&wD){
    /**
        @brief Runs the preprocessor function, if enabled.

        @tparam nodeType The type of node being processed.

        @param args The arguments to be used for preprocessor configuration.

        @param wD The pipeline packet to process.
    */
    
	//Start with the preprocessing function, if enabled
	auto pre = args["preprocessor"];
	if (pre != ""){	
		auto prePipe = preprocessor<nodeType>::newPreprocessor(pre);

		if (prePipe != 0 && prePipe->configPreprocessor(args)){
			prePipe->runPreprocessorWrapper(wD);
		}
		else{
			std::cout << "LHF processData: Failed to configure pipeline: " << args["pipeline"] << std::endl;
		}
			
		auto sv = args.find("scalarV");
		if(sv == args.end()){
			auto clusters = std::atoi(args["clusters"].c_str());
			auto scalar = std::atof(args["scalar"].c_str());
			args["scalarV"] = std::to_string(scalar * utils::computeMaxRadius(clusters, wD.workData, wD.inputData, wD.centroidLabels));
			std::cout << "Using scalarV: " << args["scalarV"] << std::endl;
		}
	}
}

template<typename nodeType>
std::vector<bettiBoundaryTableEntry> LHF<nodeType>::processParallel(std::map<std::string, std::string> args, std::vector<unsigned> &centroidLabels, std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> &partitionedData, std::vector<std::vector<double>> &inputData, int displacement){
    /**
        @brief Processes partitions in parallel and returns a merged betti table

        @tparam nodeType The node type for the pipeline

        @param args A map of arguments for pipeline configuration

        @param centroidLabels The centroid labels for each point in the data

        @param partitionedData A pair containing partition labels and partitioned data

        @param inputData The input data

        @param displacement The displacement of the first partition in the full dataset

        @return std::vector<bettiBoundaryTableEntry> The merged betti table for all partitions
    */

	//		Parameters used for each thread
	auto threshold = std::atoi(args["threshold"].c_str());
	auto maxEpsilon = std::atof(args["epsilon"].c_str());
	auto threads = std::atoi(args["threads"].c_str());
	auto clusters = std::atoi(args["clusters"].c_str());

	//		Local Storage for each thread
	std::vector<bettiBoundaryTableEntry> mergedBettiTable;
	std::vector<bettiBoundaryTableEntry> partBettiTable[threads];
	std::string runLogs[threads];
	std::string stats[threads];
	auto iterwD = pipePacket<nodeType>(args, args["complexType"]);

	//		Process fuzzy partitions in order of size
	std::vector<std::pair<unsigned, unsigned>> sortpartitions;
	for (int i = 0; i < partitionedData.second.size(); i++){
		sortpartitions.push_back(std::make_pair(partitionedData.second[i].size(), i));
	}

	//std::sort(sortpartitions.begin(), sortpartitions.end());

	std::cout << "Sorted bins: ";
	for (auto a : sortpartitions)
		std::cout << a.first << " ";

	std::cout << std::endl;

	//3. Process each partition using OpenMP to handle multithreaded scheduling
	std::cout << "Running with " << threads << " threads" << std::endl;

	//Args, but without the upscale pipe in pipeline so that we don't upscale the regional intervals
	std::map<std::string, std::string> argsNoUpscale = args;
	auto pipe = args.find("pipeline");
	if(pipe != args.end()){
		std::string pipeline = args["pipeline"];
		if(pipeline.length() > 8 && pipeline.substr(pipeline.length()-8) == ".upscale"){
			pipeline = pipeline.substr(0, pipeline.length()-8);
			argsNoUpscale["pipeline"] = pipeline;
		}
	}

	//		This begins the parallel region; # threads are spun up within these braces
	#pragma omp parallel num_threads(threads)
	{

		//		Get the thread number with OpenMP - prevents simultaneous accesses to the betti tables
		int np = omp_get_thread_num();

		//		Schedule each thread to compute one of the partitions in any order
		//			When finished, thread will grab next iteration available
		#pragma omp for schedule(dynamic)
		for (int p = partitionedData.second.size() - 1; p >= 0; p--){
			unsigned z = sortpartitions[p].second;

			//Check if we are running the centroid dataset (there's no label associated)

			//---------------------------------TODO: FIX WITH BOTH VERSIONS -------------------------------------------------------------
			if (partitionedData.second[z].size() == clusters && (args["mode"] == "mpi" || z == partitionedData.second.size() - 1)){

				std::cout << "Running centroids with " << clusters << " clusters; id = " << id << std::endl;

				//		Run on full dataset
				//		Set the pipeline to include the boundary and upscaling steps
				auto centArgs = args;
				//if(args["upscale"] == "true" || args["upscale"] == "1")
				//	centArgs["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
				//else
				//	centArgs["pipeline"] = "distMatrix.neighGraph.incrementalPersistence";

				iterwD.ident = std::to_string(np) + "," + std::to_string(p);

				//Run against the original dataset

				if (partitionedData.second[z].size() > 0){
					iterwD.inputData = inputData;
					iterwD.workData = partitionedData.second[z];
					iterwD.centroidLabels = centroidLabels;

					if(centArgs["involutedUpscale"] == "true"){
						centArgs["involuted"] = "true";
					}

					runPipeline(centArgs, iterwD);

					delete iterwD.complex;
				}
				else
					std::cout << "skipping full data, no centroids" << std::endl;

				//Determine if we need to upscale any additional boundaries based on the output of the centroid approximated PH
			} else if (sortpartitions[p].first > 0){ //Nonempty partition

				//		Clone the pipePacket<simplexNode>to prevent shared memory race conditions
				auto curwD = pipePacket<nodeType>(args, args["complexType"]);
				curwD.workData = partitionedData.second[z];
				curwD.inputData = partitionedData.second[z];
				curwD.ident = std::to_string(np) + "," + std::to_string(p);

				//		If the current partition is smaller than the threshold, process
				//			Otherwise recurse to reduce the number of points
				
				if((args["mode"] == "iter" || args["mode"] == "iterUpscale") && partitionedData.second[z].size() >= threshold){
					curwD.bettiTable = processParallelWrapper(args, curwD);
				} else{
					runPipeline(argsNoUpscale, curwD);
				}
				
				runLogs[np] += curwD.runLog;
				stats[np] += curwD.stats;

				delete curwD.complex;

				//4. Process and merge bettis - whether they are from runPipeline or IterUpscale
				unsigned h0count = 0;
				// utils::extractBoundaryPoints(curwD.bettiTable);

				//Remap the boundary indices into the original point space
				curwD.bettiTable = utils::mapPartitionIndexing(partitionedData.first[z], curwD.bettiTable);

				for (auto betEntry : curwD.bettiTable){
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
					if (betEntry.boundaryPoints.size() > 0 && centroidLabels[*boundIter] == displacement + z){
						if (betEntry.bettiDim == 0){
							boundIter++;

							//Check if second entry is in the partition
							if (centroidLabels[*boundIter] == displacement + z){
								h0count++;
								partBettiTable[np].push_back(betEntry);
							}
						}
						else{
							partBettiTable[np].push_back(betEntry);
						}
					}
				}

				unsigned additionalExternal = std::count(centroidLabels.begin(), centroidLabels.end(), displacement + z) - h0count;

				for (auto betEntry : curwD.bettiTable){
					auto boundIter = betEntry.boundaryPoints.begin();

					if (betEntry.boundaryPoints.size() > 0 && centroidLabels[*boundIter] == displacement + z){
						if (betEntry.bettiDim == 0){
							boundIter++;

							//Check if second entry is not in the partition
							if (centroidLabels[*boundIter] != displacement + z && additionalExternal > 0){
								additionalExternal--;
								partBettiTable[np].push_back(betEntry);
							}
							if (additionalExternal == 0)
								break;
						}
						else
							break;
					}
				}

				// for(auto newEntry : temp){
				// 	bool found = false;
				// 	for(auto curEntry : partBettiTable[np]){
				// 		if(newEntry.death == curEntry.death && newEntry.boundaryPoints == curEntry.boundaryPoints){
				// 			found = true;
				// 		}
				// 	}
				// 	if(!found)
				// 		partBettiTable[np].push_back(newEntry);
				// }
			}
			else std::cout << "skipping" << std::endl;
		}
	}

	//5. Merge partitions, compute PH on original centroid dataset, and report results

	//		Merge partitioned betti tables together
	for (auto partTable : partBettiTable)
		mergedBettiTable.insert(mergedBettiTable.end(), partTable.begin(), partTable.end());

	//		Merge bettis from the centroid based data
	// ---------------------------------- TODO: Fix for nprocs --------------------------------------
	if(id == nprocs - 1){
		for(auto betEntry : iterwD.bettiTable){
			betEntry.isCentroid = true;
			if(betEntry.bettiDim > 0){
				mergedBettiTable.push_back(betEntry);
			}
		}

		//		Add open d0 intervals for the remaining d0 bettis
		bettiBoundaryTableEntry des = {0, 0, maxEpsilon, {}};
		mergedBettiTable.push_back(des);
	}

	for (auto stat : stats)
		iterwD.stats += stat;

	for (auto runLog : runLogs)
		iterwD.runLog += runLog;

	iterwD.bettiTable = mergedBettiTable;
	outputBettis(args, iterwD);

	//		Return the final merged betti table for this iteration
	return mergedBettiTable;
}


template<typename nodeType>
std::vector<bettiBoundaryTableEntry> LHF<nodeType>::processParallelWrapper(std::map<std::string, std::string> args, pipePacket<nodeType>&wD, bool runPartition){

	//This function is called when the number of points in a partition are greater than the point threshold
	//	If the number of points in the new partitions are under the point threshold, continue
	//	If the number of points in any new partition is above the point threshold, recurse
	//
	//We may also want to limit the number of recursions for bigger data sets

	//Break this down into better code sections:

	//1. Local objects and storage

	//		Parameters
	auto scalar = std::atof(args["scalar"].c_str());
	auto clusters = std::atoi(args["clusters"].c_str());

	//		Local Storage
	auto originalDataSize = wD.inputData.size();

	//2. Partition the source point cloud separate datasets accordingly
	if (runPartition){
		//		Run preprocessor pipeline to partition
		runPreprocessor(args, wD);
	}

	//		Compute partition statistics for fuzzy partition distance
	auto maxRadius = utils::computeMaxRadius(clusters, wD.workData, wD.inputData, wD.centroidLabels);
	//auto avgRadius = utils::computeAvgRadius(clusters, wD.workData, wD.inputData, wD.centroidLabels);

	//		Count the size of each partition for identifying source partitions when looking at the betti table results
	std::vector<unsigned> binCounts;
	for (unsigned a = 0; a < clusters; a++){
		binCounts.push_back(std::count(wD.centroidLabels.begin(), wD.centroidLabels.end(), a));
	}
	std::cout << "Bin Counts: ";
	utils::print1DVector(binCounts);

	//		Sort our fuzzy partitions into individual vectors
	args["scalarV"] = std::to_string(scalar * maxRadius);
	auto partitionedData = utils::separatePartitions(std::atof(args["scalarV"].c_str()), wD.workData, wD.inputData, wD.centroidLabels);
	std::cout << "Using scalar value: " << args["scalarV"] << std::endl;

	//		Append the centroid dataset to run in parallel as well
	partitionedData.second.push_back(wD.workData);

	return processParallel(args, wD.centroidLabels, partitionedData, wD.inputData);
}


template<typename nodeType>
std::vector<bettiBoundaryTableEntry> LHF<nodeType>::processDistributedWrapper(std::map<std::string, std::string> args, pipePacket<nodeType>&wD){
	
	//Local arguments for controlling partitioning and merging
	auto scalar = std::atof(args["scalar"].c_str());
	auto clusters = std::atoi(args["clusters"].c_str());

	//Local classes for reading, writing, utilities
	auto rs = readInput();

	// check if 1 process only
	if (nprocs == 1){
		//Read our input data
		wD.inputData = rs.readCSV(args["inputFile"]);
		wD.workData = wD.inputData;

		return processParallelWrapper(args,wD);
	}

	//Final Betties result will be saved here
	std::vector<bettiBoundaryTableEntry> finalMergedBettiTable;

	//local Buffers to serialize betties in different processes
	std::vector<unsigned> betti_dim;
	std::vector<double> betti_birth;
	std::vector<double> betti_death;
	std::vector<unsigned> betti_boundarysize;
	std::vector<unsigned> betti_boundaries;

	//Data dimension
	unsigned dimension;

	//Minimum Partitions per process
	int minPartitions = 0;

	//During distribution first k processes that will receive one more partition than minPartitions
	int firstk = 0;

	//Each Process is aware of fuzzy partition sizes
	std::vector<unsigned> partitionsize;
	unsigned partitionsize_size = 0;

	//Serialised Point Cloud and their labels
	std::vector<double> sdata;
	std::vector<unsigned> sdatalabel;

	//Each process is aware of original labels
	std::vector<unsigned> originalLabels;
	unsigned originalLabels_size;

	//Chunk size and offsets of data and their labels for MPI functions
	int scountslabel[nprocs];
	int displslabel[nprocs];
	int scounts[nprocs];
	int displs[nprocs];

	//Buffers and their sizes for receiving chunks
	int maxsize = 0;
	int maxsizelabel = 0;
	double *receivedData;
	unsigned *receivedDataLabel;

	//Betties table sizes and their boundary sizes generated by each processes
	int bettiTableSize = 0;
	unsigned boundary_size = 0;

	//Check if we are the master process for upscaling
	if (id == 0){

		//Read our input data
		wD.inputData = rs.readCSV(args["inputFile"]);
		wD.workData = wD.inputData;

		//Partition the data with the configured preprocessor
		runPreprocessor(args, wD);

		//original Labels to be sent our to all processes
		originalLabels = wD.centroidLabels;
		originalLabels_size = wD.centroidLabels.size();

		//Separate our partitions for distribution
		auto maxRadius = utils::computeMaxRadius(std::atoi(args["clusters"].c_str()), wD.workData, wD.inputData, wD.centroidLabels);
		auto avgRadius = utils::computeAvgRadius(std::atoi(args["clusters"].c_str()), wD.workData, wD.inputData, wD.centroidLabels);

		auto partitionedData1 = utils::separatePartitions(scalar * maxRadius, wD.workData, wD.inputData, wD.centroidLabels);

		partitionedData1.second.push_back(wD.workData);
		//	Each node/slave will process at least 1 partition
		//	NOTE: the partition may contain points outside partition that are within 2*Rmax
		minPartitions = partitionedData1.second.size() / nprocs;
		firstk = partitionedData1.second.size() - (minPartitions * nprocs);

		dimension = partitionedData1.second[0][0].size();
		//Get the partition sizes
		for (auto a : partitionedData1.second)
			partitionsize.push_back(a.size());

		//serializing all the data
		for (auto a : partitionedData1.second)
			for (auto b : a)
				for (auto c : b)
					sdata.push_back(c);

		//serialize all the labels
		for (auto a : partitionedData1.first)
			for (auto b : a)
				sdatalabel.push_back(b);

		//intitalizing to 0
		for (int i = 0; i < nprocs; i++){
			scounts[i] = 0;
			displs[i] = 0;
			scountslabel[i] = 0;
			displslabel[i] = 0;
		}

		//distributing partitions sizes in scounts and scountslabel
		int temp = 0, k = 0;
		for (auto a : partitionsize){
			if (k < firstk){
				if (temp >= minPartitions + 1){
					temp = 0;
					k++;
				}
			}
			else{
				if (temp >= minPartitions){
					temp = 0;
					k++;
				}
			}
			scounts[k] += a * dimension;
			scountslabel[k] += a;
			temp++;
		}

		k = 0;
		maxsize = 0;
		maxsizelabel = 0;
		auto bd = scounts[0];

		// find the maximum buffer size required for partition + offset for each partition
		for (auto a : scounts){
			if (k > 0)
				displs[k] = bd + displs[k - 1];
			k++;
			bd = a;

			if (a > maxsize)
				maxsize = a;
		}

		k = 0;
		auto bl = scountslabel[0];
		// find maximum buffer size required by label partition	+ offset for each partition
		for (auto a : scountslabel){
			if (k > 0)
				displslabel[k] = bl + displslabel[k - 1];
			k++;
			bl = a;

			if (a > maxsizelabel)
				maxsizelabel = a;
		}

		partitionsize_size = partitionsize.size();
		std::cout << "partitionsize:";
		utils::print1DVector(partitionsize);
	}

	//	broadcasting all the required information by slaves.

	MPI_Bcast(&dimension, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);			//Data dimension
	MPI_Bcast(&partitionsize_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD); //Total number of partitions
	partitionsize.resize(partitionsize_size);
	MPI_Bcast(&partitionsize[0], partitionsize_size, MPI_UNSIGNED, 0, MPI_COMM_WORLD); //each partition size array

	MPI_Bcast(&originalLabels_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD); //total number of labels
	originalLabels.resize(originalLabels_size);
	MPI_Bcast(&originalLabels[0], originalLabels_size, MPI_UNSIGNED, 0, MPI_COMM_WORLD); //label array

	//Partitions specific details
	MPI_Bcast(&firstk, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&minPartitions, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Maximum buffer required information
	MPI_Bcast(&maxsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxsizelabel, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Chunk sizes for each processes
	MPI_Bcast(&scounts, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&scountslabel, nprocs, MPI_INT, 0, MPI_COMM_WORLD);

	//firstk receive one more partition
	if (id < firstk)
		minPartitions = minPartitions + 1;

	// Allocating buffers
	receivedData = (double *)malloc(maxsize * sizeof(double));
	receivedDataLabel = (unsigned *)malloc(maxsizelabel * sizeof(unsigned));
	//Initializing buffers to zeros
	for (int i = 0; i < maxsize; i++){
		receivedData[i] = 0.0;
	}
	for (int i = 0; i < maxsizelabel; i++){
		receivedDataLabel[i] = 0;
	}

	//scatter the point cloud partitions to slaves.
	MPI_Scatterv(&sdata[0], scounts, displs, MPI_DOUBLE, receivedData, maxsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// scatter label information
	MPI_Scatterv(&sdatalabel[0], scountslabel, displslabel, MPI_UNSIGNED, receivedDataLabel, maxsizelabel, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if (minPartitions > 0){

		int displacement = 0; // offset for each processes to read partitions
		if (id < firstk)
			displacement = id * minPartitions;
		else
			displacement = firstk * (minPartitions + 1) + (id - firstk) * minPartitions;

		//each process separates their partitions and labels
		std::vector<std::vector<std::vector<double>>> partitionData;
		std::vector<std::vector<unsigned>> labels;

		int counter = 0; //counter

		//Each process separate their data partitions

		for (int i = 0; i < minPartitions; i++){
			std::vector<std::vector<double>> partition;
			for (int j = 0; j < partitionsize[displacement + i]; j++){
				std::vector<double> row;
				for (int k = 0; k < dimension; k++){
					row.push_back(receivedData[counter++]);
				}

				partition.push_back(row);
			}

			partitionData.push_back(partition);
		}

		//Each process separate thier label partitions
		counter = 0;
		for (int i = 0; i < minPartitions; i++){
			std::vector<unsigned> partitionlabel;
			for (int j = 0; j < partitionsize[displacement + i]; j++){
				partitionlabel.push_back(receivedDataLabel[counter++]);
			}
			labels.push_back(partitionlabel);
		}

		auto partitionedData = make_pair(labels, partitionData);

		//Local Storage
		auto originalDataSize = wD.workData.size();

		std::vector<bettiBoundaryTableEntry> mergedBettiTable = processParallel(args, originalLabels, partitionedData, wD.inputData, displacement);

		//Serialize bettie table to send to master process for merging
		for (auto bet : mergedBettiTable){
			betti_dim.push_back(bet.bettiDim);
			betti_birth.push_back(bet.birth);
			betti_death.push_back(bet.death);
			betti_boundarysize.push_back(bet.boundaryPoints.size());
			betti_boundaries.insert(betti_boundaries.end(), bet.boundaryPoints.begin(), bet.boundaryPoints.end());
		}

		bettiTableSize = betti_dim.size();
		boundary_size = 0;
		for (auto a : betti_boundarysize)
			boundary_size += a;
	}
	else{
		std::cout << "**********************\n"
				  << std::endl;
	}

	//gather all the betties from slaves
	int totalsize = 0;
	//Sum up all the betties sizes from diffrent processes at process 0
	MPI_Reduce(&bettiTableSize, &totalsize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	//Buffers in which all the betties for diffrent processes will be collected
	std::vector<unsigned> recvbuffdim(totalsize);
	std::vector<double> recvbuffbirth(totalsize);
	std::vector<double> recvbuffdeath(totalsize);
	std::vector<unsigned> recvbuffboundaries_size(totalsize);

	//gether each process table sizes
	std::vector<int> betties_table_size(nprocs);
	MPI_Gather(&bettiTableSize, 1, MPI_INT, &betties_table_size[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

	std::vector<int> displsg(nprocs);
	displsg[0] = 0;

	// Calculate offsets
	for (int i = 1; i <= nprocs; i++)
		displsg[i] = displsg[i - 1] + betties_table_size[i - 1];

	// Gather betties info
	MPI_Gatherv(&betti_dim[0], bettiTableSize, MPI_UNSIGNED, &recvbuffdim[0], &betties_table_size[0], &displsg[0], MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Gatherv(&betti_birth[0], bettiTableSize, MPI_DOUBLE, &recvbuffbirth[0], &betties_table_size[0], &displsg[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(&betti_death[0], bettiTableSize, MPI_DOUBLE, &recvbuffdeath[0], &betties_table_size[0], &displsg[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(&betti_boundarysize[0], bettiTableSize, MPI_UNSIGNED, &recvbuffboundaries_size[0], &betties_table_size[0], &displsg[0], MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	int totalboundarysize = 0;
	MPI_Reduce(&boundary_size, &totalboundarysize, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

	//Gather boundary size information
	std::vector<int> betties_table_boundary_size(nprocs);
	MPI_Gather(&boundary_size, 1, MPI_UNSIGNED, &betties_table_boundary_size[0], 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	std::vector<int> displsgb(nprocs);
	displsgb[0] = 0;

	//Calculate offsets
	for (int i = 1; i <= nprocs; i++){
		displsgb[i] = displsgb[i - 1] + betties_table_boundary_size[i - 1];
	}

	std::vector<unsigned> recvbuffboundaries(totalboundarysize);
	//Gather Boundary information
	MPI_Gatherv(&betti_boundaries[0], boundary_size, MPI_UNSIGNED, &recvbuffboundaries[0], &betties_table_boundary_size[0], &displsgb[0], MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if (id == 0){
		//master prune out the duplicates across slaves.

		int beg = 0;
		for (int i = 0; i < totalsize; i++){
			bettiBoundaryTableEntry bettiEntry;
			bettiEntry.bettiDim = recvbuffdim[i];
			bettiEntry.birth = recvbuffbirth[i];
			bettiEntry.death = recvbuffdeath[i];
			for (int bi = beg; bi < (beg + recvbuffboundaries_size[i]); bi++)
				bettiEntry.boundaryPoints.insert(recvbuffboundaries[bi]);

			beg += recvbuffboundaries_size[i];
			bool found = false;
			for (auto curEntry : finalMergedBettiTable){
				if (bettiEntry.death == curEntry.death && bettiEntry.boundaryPoints == curEntry.boundaryPoints){
					found = true;
				}
			}

			if (!found)
				finalMergedBettiTable.push_back(bettiEntry);
		}
	}

	return finalMergedBettiTable;
}

template<typename nodeType>
pipePacket<nodeType> *processpyLHFWrapper(std::map<std::string, std::string> &args, std::vector<std::vector<double>> &pointCloud){
    
    //Determine what pipe we will be running; this involves creating a temporary LHF
    //  object and using a void pointer to reassign the correct complex type
    auto lhflib = LHF<nodeType>();
    
    auto wD = new pipePacket<nodeType>(args, args["complexType"]);

    wD->inputData = pointCloud;
    wD->workData = wD->inputData;
    
    double start = omp_get_wtime();
    
    if (wD->inputData.size() > 0 || args["pipeline"] == "slidingwindow" || args["pipeline"] == "naivewindow" || args["mode"] == "mpi"){

        if(args["mode"] == "reduced" || args["mode"] == "iterUpscale" || args["mode"] == "iter"){	
            wD->bettiTable = lhflib.processParallelWrapper(args,*wD);
            sort(wD->bettiTable.begin(), wD->bettiTable.end(), sortBettis());
        }
        else{
            lhflib.runPreprocessor(args, *wD);
            lhflib.runPipeline(args, *wD);
        }
    }
    else{
        argParser::printUsage();
    }

    if ((args["debug"] == "1" || args["debug"] == "true") && wD->bettiTable.size() > 0){
        std::cout << std::endl
                  << "_______Merged BETTIS_______" << std::endl;

        for (auto a : wD->bettiTable){
            std::cout << a.bettiDim << ",\t" << a.birth << ",\t" << a.death << ",\t";
            utils::print1DVector(a.boundaryPoints);
        }
    }

    delete wD->complex;

    double end = omp_get_wtime();
    std::cout << "Total LHF execution time (s): " << end - start << std::endl;
    
    return wD;
}



extern "C"{
	
    
	pipeWrap *pyLHFWrapper(int argc, char *argv, const double *pointCloud){
        pipeWrap *ret;
        
        /*******     1. Read Arguments      *********/
		std::map<std::string, std::string> args;
		std::vector<std::string> rawArgs;

		//Split arguments into list
		std::string tempstr = "";
		for (auto i = 0; i < argc; i++){
			//Check for space (32)
			if (argv[i] == 32){
				rawArgs.push_back(tempstr);
				tempstr = "";
			}
			else
				tempstr += argv[i];
		}

		//Split argument list into map
		for (auto i = 0; i < rawArgs.size(); i += 2){
			args[rawArgs[i]] = rawArgs[i + 1];
		}
        
        
        
        //Get args from pipeline mode
        argParser::defaultArguments(args);
        argParser::setPipeline(args);

		/*******     2. Decode Data          *********/
		int dataSize = std::atoi(args["datasize"].c_str());
		int dataDim = std::atoi(args["datadim"].c_str());

		std::vector<std::vector<double>> data(dataSize, std::vector<double>(dataDim));

		for (auto row = 0; row < dataSize; row++){
			for (auto dim = 0; dim < dataDim; dim++){
				data[row][dim] = pointCloud[row*dataDim + dim];
			}
		}
        
        
		/*******     3. Run LHF             *********/
        //Local pointers for different nodeTypes
        std::vector<bettiBoundaryTableEntry>* l_bettiTable;
        std::string *l_ident;
        std::string *l_stats;
        std::string *l_runLog;
        std::vector<std::vector<double>> *l_workData;
        std::vector<unsigned> *l_centroidLabels;
        std::vector<std::vector<double>> *l_inputData;
        std::vector<std::vector<double>> *l_distMatrix;
        std::vector<std::vector<bool>> *l_incidenceMatrix;
        std::vector<std::set<unsigned>> *l_boundaries;
        
        if(args["nodeType"] == "alphaNode"){
            auto ret = processpyLHFWrapper<alphaNode>(args, data);
            l_bettiTable = &ret->bettiTable;
            l_ident = &ret->ident;
            l_stats = &ret->stats;
            l_runLog = &ret->runLog;
            l_workData = &ret->workData;
            l_centroidLabels = &ret->centroidLabels;
            l_inputData = &ret->inputData;
            l_distMatrix = &ret->distMatrix;
            l_incidenceMatrix = &ret->incidenceMatrix;
            l_boundaries = &ret->boundaries;
            
        } else if (args["nodeType"] == "witnessNode"){
            auto ret = processpyLHFWrapper<witnessNode>(args, data);
            l_bettiTable = &ret->bettiTable;
            l_ident = &ret->ident;
            l_stats = &ret->stats;
            l_runLog = &ret->runLog;
            l_workData = &ret->workData;
            l_centroidLabels = &ret->centroidLabels;
            l_inputData = &ret->inputData;
            l_distMatrix = &ret->distMatrix;
            l_incidenceMatrix = &ret->incidenceMatrix;
            l_boundaries = &ret->boundaries;
        } else {
            auto ret = processpyLHFWrapper<simplexNode>(args, data);
            l_bettiTable = &ret->bettiTable;
            l_ident = &ret->ident;
            l_stats = &ret->stats;
            l_runLog = &ret->runLog;
            l_workData = &ret->workData;
            l_centroidLabels = &ret->centroidLabels;
            l_inputData = &ret->inputData;
            l_distMatrix = &ret->distMatrix;
            l_incidenceMatrix = &ret->incidenceMatrix;
            l_boundaries = &ret->boundaries;
        }
        
        
        retPipePacket* retStruct = (retPipePacket *) malloc (sizeof(retPipePacket));
        
        //Calculate the size of the betti table with total generator entries
        int generators = 0;
        for(auto i = 0; i < (*l_bettiTable).size(); i++){
            generators += (*l_bettiTable)[i].boundaryPoints.size();
        }
        
        //Allocate and populate the serialized betti table; size is number of entries plus the total number of generators
		retBettiTable *bettiTable = (retBettiTable *)malloc(sizeof(retBettiTable) * (*l_bettiTable).size() + sizeof(unsigned) * generators); 

		for (auto i = 0; i < (*l_bettiTable).size(); i++){
			bettiTable[i].dim = (*l_bettiTable)[i].bettiDim;
			bettiTable[i].birth = (*l_bettiTable)[i].birth;
			bettiTable[i].death = (*l_bettiTable)[i].death;
            bettiTable[i].boundarySize = (*l_bettiTable)[i].boundaryPoints.size();
            
            bettiTable[i].boundaryEntries = (unsigned*)malloc(sizeof(unsigned) * bettiTable[i].boundarySize);
            
            int jdx = 0;
            for (auto generator: (*l_bettiTable)[i].boundaryPoints){
                bettiTable[i].boundaryEntries[jdx] = generator;
                jdx++;
            }
            
		}
        
        retStruct->size_betti = (*l_bettiTable).size();
        retStruct->LHF_size = dataSize;
        retStruct->LHF_dim = dataDim;
        retStruct->workData_size = (*l_workData).size();
        retStruct->bettiTable = bettiTable;                 
        retStruct->stats = (*l_stats).c_str();               
        retStruct->runLog = (*l_runLog).c_str();
        retStruct->ident = (*l_ident).c_str();
        
        double *inputData_retStruct = (double *)malloc(sizeof(double) * ((*l_inputData).size() * (*l_inputData)[0].size()));
        double *distMatrix_retStruct = (double *)malloc(sizeof(double) * ((*l_distMatrix).size() * (*l_distMatrix).size()));
        double *workData_retStruct = (double *)malloc(sizeof(double) * ((*l_workData).size() * (*l_workData).size()));
        unsigned *centroidLabels_retStruct = (unsigned *)malloc(sizeof(unsigned) * ((*l_centroidLabels).size()));
		//char *stats_retStruct = (char *)malloc(sizeof(char) * wD.stats.length());
		//char *runLog_retStruct = (char *)malloc(sizeof(char) * wD.runLog.length());
		//char *ident_retStruct = (char *)malloc(sizeof(char) * wD.ident.length());
        
		auto wD = pipePacket<simplexNode>(args, args["complexType"]); //wD (workingData)
        
        inputData_retStruct = &std::accumulate((*l_inputData).begin(), (*l_inputData).end(), decltype(wD.inputData)::value_type{}, [](auto &x, auto &y){x.insert(x.end(), y.begin(), y.end()); return x;})[0];
        distMatrix_retStruct = &std::accumulate((*l_distMatrix).begin(), (*l_distMatrix).end(), decltype(wD.distMatrix)::value_type{}, [](auto &x, auto &y){x.insert(x.end(), y.begin(), y.end()); return x;})[0];
        workData_retStruct = &std::accumulate((*l_workData).begin(), (*l_workData).end(), decltype(wD.workData)::value_type{}, [](auto &x, auto &y){x.insert(x.end(), y.begin(), y.end()); return x;})[0];
        centroidLabels_retStruct = &((*l_centroidLabels).begin())[0];
        
        retStruct->inputData = inputData_retStruct;
        retStruct->distMatrix = distMatrix_retStruct;
        retStruct->centroidLabels = centroidLabels_retStruct;
        retStruct->workData = workData_retStruct;

		std::cout << "Finished on the c++ side; Pis -> " << (*l_bettiTable).size() << std::endl;

        
        return retStruct;
    }


}

