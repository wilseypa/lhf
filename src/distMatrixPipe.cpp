/*
 * distMatrix hpp + cpp extend the basePipe class for calculating the 
 * distance matrix from data input
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include "distMatrixPipe.hpp"

// basePipe constructor
distMatrixPipe::distMatrixPipe(){

	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket distMatrixPipe::runPipe(pipePacket inData){
		
	std::cout << "Made it to distMatrixPipe::runPipe" << std::endl;
	
	//inDat = std::sqrt(std::inner_product(inData.begin(), inData.end(), enddata.begin,
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool distMatrixPipe::configPipe(std::map<std::string, std::string> configMap){
	
	std::cout << "Made it to distMatrixPipe::configPipe" << std::endl;
	return true;
}

