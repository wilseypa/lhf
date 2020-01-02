/*
 * writeOutput hpp + cpp protoype and define a class for writing output from
 * the LFH system (https://github.com/wilseypa/LFH). The class is designed
 * to hold all functions corresponding to wiriting output data to .csv files,
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

bool writeOutput::writeBarcodes(std::vector<std::vector<double>> data, std::string filename){
	std::string header = "dimension,birth,death\n";
	
	return writeCSV(data, filename, header);
}

bool writeOutput::writeBarcodes(std::string data, std::string filename){
	std::string header = "dimension,birth,death\n";
	
	return writeCSV(data, filename, header);
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
bool writeOutput::writeMAT(std::vector<std::vector<double>>, std::string filename){
	std::vector<std::vector<double>> result;
	
	std::ifstream file;
	file.open(filename);
	
	std::string line;			//Temporary (current) line
	
	// Get the number of vectors
	getline(file,line);
	line = std::regex_replace(line,std::regex(" "),"");
	int vectors = std::stoi(line);
	
	// Get the number of dimensions
	getline(file,line);
	line = std::regex_replace(line,std::regex(" "),"");
	int dimensions = std::stoi(line);
	
	// We are going to iterate through each line of the file until we reach the end
	while(!file.eof()){
		std::vector<double> tmp;	//Temporary (current) vector
		getline(file, line);		//Read the next line from file
		std::size_t p = std::string::npos;
		
		// Replace whitespace in the current line
		line = std::regex_replace(line, std::regex(" "), "");
		
		// Check if the line has a length (is not a blank line)
		if(line.size() > 0){
			
			// Iterate through each comma of the csv
			while((p = line.find_first_of(",")) != std::string::npos){
				
				// Push the value found before the comma, remove from the line
				tmp.push_back(std::stod(line.substr(0,p)));
				line.erase(0,p + 1);
				
			}
			
			// Get the last value of the line
			tmp.push_back(std::stod(line.substr(0,p)));
		}
		result.push_back(tmp);
		
	}
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
