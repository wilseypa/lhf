/*
 * readInput hpp + cpp protoype and define a class for reading input into
 * the LFH system (https://github.com/wilseypa/LFH). The class is designed
 * to hold all functions corresponding to reading input data from .csv files,
 * .mat files, and other additional files as they are added to the system.
 * 
 *  
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include <regex>
#include <fstream>
#include "readInput.hpp"

// readInput constructor, currently no needed information for the class constructor
readInput::readInput(){

}

// readCSV -> read in a csv formatted file from file input
//		-filename - complete filename and relative path (if needed) for reading
// 
// Formatted as: 1.0423, 1.0244, 1.032 \n
//
// By default, double conversion (std::stod) will handle 'E' and 'e' for scientific notation
// 
// Strips whitespace characters where appropriate to create a vector array of doubles
//
// TODO: Check for invalid vector lengths (i.e. <3,3,3> and <4,4,4,4> should not exist)
// TODO: Handle if there is a comma at the end of the line (i.e. "5,5,5," should be <5,5,5>)
std::vector<std::vector<double>> readInput::readCSV(std::string filename){
	std::vector<std::vector<double>> result;
	
	std::ifstream file;
	file.open(filename);
	
	// We are going to iterate through each line of the file until we reach the end
	while(!file.eof()){
		std::string line;			// Temporary (current) line
		std::vector<double> tmp;	// Temporary (current) vector
		getline(file, line);		// Read the next line from file
		std::size_t pos = std::string::npos;
		
		// Replace whitespace in the current line
		line = std::regex_replace(line, std::regex(" "), "");
		
		// Check if the line has a length (is not a blank line)
		if(line.size() > 0){
			
			// Iterate through each comma of the csv
			while((pos = line.find_first_of(",")) != std::string::npos){
				
				// Push the value found before the comma, remove from the line
				tmp.push_back(std::stod(line.substr(0,pos)));
				line.erase(0,pos + 1);
				
			}
			
			// Get the last value of the line
			tmp.push_back(std::stod(line.substr(0,pos)));
		}
		result.push_back(tmp);
		
	}
	return result;
}


// readMAT -> read in a mat formatted file from file input
//		-filename - complete filename and relative path (if needed) for reading
// 
// Formatted as: 	15		(# of vectors)
//					20		(# of dimensions)
//					1.023	(value of [0,0])
//					4.234	(value of [0,1])
//					...		(subsequent values of [i, j])
//
// By default, double conversion (std::stod) will handle 'E' and 'e' for scientific notation
// 
// Strips whitespace characters where appropriate to create a vector array of doubles
//
// TODO: A lot
std::vector<std::vector<double>> readInput::readMAT(std::string filename){
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
	return result;
}

