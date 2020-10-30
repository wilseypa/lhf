/*
 * writeOutput hpp + cpp protoype and define a class for writing output from
 * the LHF system (https://github.com/wilseypa/LHF). The class is designed
 * to hold all functions corresponding to writing output data to .csv files,
 * .mat files, console, and other additional files as they are added to the system.
 * 
 *  
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include <regex>
#include <fstream>
#include "writeOutput.hpp"

// writeOutput constructor, currently no needed information for the class constructor
writeOutput::writeOutput(){

}

// writeStat -> write the pipeline statistics to a csv formatted file
bool writeOutput::writeStats(std::string stats, std::string filename){
	std::ofstream file(filename + "_stats.csv");
	file << "PipeName,PhysExecutionTime,MemorySize,MemoryUnits,VertexCount,SimplexCount\n";
	file << stats;
	file.close();
	return true;	
}

bool writeOutput::writeBarcodes(std::vector<bettiBoundaryTableEntry> data, std::string filename){
	std::string header = "dimension,birth,death\n";
	std::ofstream file(filename + ".csv");
	
	for(auto row : data)
		file << std::to_string(row.bettiDim) << "," << std::to_string(row.birth) << "," << std::to_string(row.death) << std::endl;
	
	file.close();
	return true;
}

// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(std::vector<std::vector<double>> data, std::string filename){
	std::ofstream file(filename + ".csv");
	
	for(auto row : data){
		
		for(auto column : row){
			file << std::to_string(column) << ",";
		}
		file << "\n";
	}
	file.close();
	return true;
}

// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(std::string data, std::string filename){
	std::ofstream file(filename + ".csv");
	file << data;
	file.close();
	return true;
}

// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(std::vector<std::vector<double>> data, std::string filename, std::string header){
	std::ofstream file(filename + ".csv");
	
	file << header;
	
	for(auto row : data){
		
		for(auto column : row){
			file << std::to_string(column) << ",";
		}
		file << "\n";
	}
	file.close();
	return true;
	
}

// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(std::string data, std::string filename, std::string header){
	std::ofstream file(filename + ".csv");
	file << header;
	file << data;
	file.close();
	return true;
	
}


// writeMAT -> write in a mat formatted file of data input
bool writeOutput::writeMAT(std::vector<std::vector<double>> data, std::string filename){
	std::ofstream file(filename + ".mat");
	
	file << std::to_string(data.size()) << "\n" << std::to_string(data[0].size()) << "\n";
	for(auto row : data){
		
		for(auto column : row){
			file << std::to_string(column) << "\n";
		}
	}
	file.close();
	return true;
}

/*
// writeConsole -> write data input to console (hopefully pretty-print)
bool writeOutput::writeConsole(pipePacket* workData){
	std::vector<std::vector<double>> result;
	
	std::cout << "_________Persistent Homology OUTPUT__________" << std::endl;
		
	std::cout << workData->bettiOutput << std::endl;
	
	return true;
	
}
*/
