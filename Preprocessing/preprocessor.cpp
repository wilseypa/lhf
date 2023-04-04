/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute  
 * 
 */

/**
	@file
	@brief This file includes the headers for preprocessor.hpp and kMeansPlusPlus.hpp.
*/
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "preprocessor.hpp"
#include "kMeansPlusPlus.hpp"

/*
#include "streamingKmeans.hpp"
#include "denStream.hpp"
*/


// basePipe constructor
/**
	@brief Default constructor for the preprocessor class.
	This constructor initializes the preprocessor object.
*/
template<typename nodeType>
preprocessor<nodeType>::preprocessor(){
	return;
}

/**
	@brief Constructor for preprocessor class.
	@tparam nodeType The type of data stored in nodes of the tree.
*/

template<typename nodeType>
preprocessor<nodeType>* preprocessor<nodeType>::newPreprocessor(const std::string &procName){
	utils ut;
	ut.writeDebug("preprocessor","Building preprocessor: " + procName);

/**
	@brief Returns a pointer to a new preprocessor object based on the specified name.
	@tparam nodeType The type of data stored in nodes of the tree.
	@param procName A string indicating the type of preprocessor to create.
	@return A pointer to a new preprocessor object.
*/

	if(procName == "none"){
		return new preprocessor<nodeType>();
	} else if (procName == "kmeansplusplus" || procName == "kmeans++" || procName == "kmeans"){
		return new kMeansPlusPlus<nodeType>();
	} 
	
	/*else if(procName == "streamingKmeans" || procName == "streamingkmeans" || procName =="streamKM"){
		return new streamingKmeans();
	} else if(procName == "denStream" || procName == "denstream" || procName =="DenStream"){
		return new denStream();
	} */

	return 0;
}

/**

	@brief Wrapper function for timing of runPreprocessor and other misc. functions.

	@tparam nodeType The type of nodes in the data set.

	@param inData The input data packet for the preprocessor.

	This function checks if the preprocessor has been configured and then either executes the runPreprocessor function directly

	or starts a timer for physical time passed during the runPreprocessor function. The timer is stopped after the function has

	completed and the duration is calculated. The function then outputs the time and memory used for this pipeline segment.

	@note If debug mode is enabled, the output includes time and memory usage information.
*/
// runPipeWrapper -> wrapper for timing of runPipe and other misc. functions
template<typename nodeType>
void preprocessor<nodeType>::runPreprocessorWrapper(pipePacket<nodeType> &inData){
	
	//Check if the preprocessor has been configured
	if(!configured){
		ut.writeLog(procName,"Pipe not configured");
		std::cout << "Pipe not configured" << std::endl;
		return;
	}
	
	if(debug){
		
		//Start a timer for physical time passed during the pipe's function
		auto startTime = std::chrono::high_resolution_clock::now();
		
		runPreprocessor(inData);
		
		//Stop the timer for time passed during the pipe's function
		auto endTime = std::chrono::high_resolution_clock::now();
		
		//Calculate the duration (physical time) for the pipe's function
		std::chrono::duration<double, std::milli> elapsed = endTime - startTime;
		
		//Output the time and memory used for this pipeline segment
		std::cout << "\tPipeline " << procName << " executed in " << (elapsed.count()/1000.0) << " seconds (physical time)" << std::endl << std::endl;
		
		auto dataSize = inData.getSize();
		auto unit = "B";
		
		std::cout << "Test" << std::endl;
		//Convert large datatypes (GB, MB, KB)
		if(dataSize > 1000000000){
			//Convert to GB
			dataSize = dataSize/1000000000;
			unit = "GB";
		} else if(dataSize > 1000000){
			//Convert to MB
			dataSize = dataSize/1000000;
			unit = "MB";
		} else if (dataSize > 1000){
			//Convert to KB
			dataSize = dataSize/1000;
			unit = "KB";
		}
		
		inData.stats += procName + "," + std::to_string(elapsed.count()/1000.0) + "," + std::to_string(dataSize) + "," + unit + "\n";
		
		outputData(inData);
	
	} else {
		runPreprocessor(inData);
	}
}
/**

	@brief Outputs the given data to a CSV file in the output directory

	@tparam nodeType the type of node being processed

	@param data a vector of unsigned integers to output to the CSV file

	@return void
*/

template<typename nodeType>
void preprocessor<nodeType>::outputData(std::vector<unsigned> data){
	std::ofstream file;
	file.open("output/" + procName + "_label_output.csv");
	
	for (auto a : data){
		file << std::to_string(a) << "\n";
	}
	file.close();
	return;
}
	
/**

	@brief Outputs data to file.

	@tparam nodeType Data type of nodes.

	@param data The data to output.
*/

template<typename nodeType>
void preprocessor<nodeType>::outputData(pipePacket<nodeType> &data){
	
	outputData(data.workData);
	outputData(data.centroidLabels);
	
	return;
	
	
}

/**

	@brief Outputs a vector of vectors to a CSV file in the output folder with the name of the preprocessor and "_centroid_output.csv" appended

	@tparam nodeType the type of the nodes in the vector

	@param data the vector of vectors to output to CSV

	@return void
*/
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template<typename nodeType>
void preprocessor<nodeType>::outputData(std::vector<std::vector<double>> data){
	std::ofstream file;
	file.open("output/" + procName + "_centroid_output.csv");
	
	for (auto a : data){
		for (auto d : a){
			file << std::to_string(d) << ",";
		}
		file << "\n";
	}
	
	file.close();
	return;
}
	
// runPipe -> Run the configured functions of this pipeline segment

/**

	@brief A default implementation of the runPreprocessor function that displays an error message if called.

	@tparam nodeType The datatype of the input data nodes.

	@param inData A reference to the input data packet.
*/
template<typename nodeType>
void preprocessor<nodeType>::runPreprocessor(pipePacket<nodeType> &inData){
	
	std::cout << "No run function defined for: " << procName << std::endl;
	
	return;
}	

// configPipe -> configure the function settings of this pipeline segment
/**
	@brief Configures the preprocessor with the specified configuration options.
	@tparam nodeType The datatype of the nodes in the data structure being processed.
	@param configMap A map of string keys and values representing the configuration options.
	@return true if the preprocessor is successfully configured, false otherwise.
*/
template<typename nodeType>
bool preprocessor<nodeType>::configPreprocessor(std::map<std::string, std::string> &configMap){
	std::cout << "No configure function defined for: " << procName << std::endl;
	
	return false;
}

/**

	@brief Preprocessor class implementation for preprocessing data in the TDA pipeline.

	This template class is used to preprocess data in the TDA pipeline by running a variety of

	preprocessor functions on the input data. The template class is parameterized by nodeType,

	which determines the type of node used to represent data in the pipeline.

	@tparam nodeType Type of node used to represent data in the pipeline.
*/

template<typename nodeType>
preprocessor<nodeType>::~preprocessor(){}

//Explicit Template Class Instantiation
template class preprocessor<simplexNode>;
template class preprocessor<alphaNode>;
template class preprocessor<witnessNode>;
