#include <iostream>
#include <stdio.h>
#include <vector>
#include "writeOutput.hpp"
#include "readInput.hpp"
#include "utils.hpp"

// TEST write Functions
void t_write_functions(std::string &log){
	auto *ws = new writeOutput();
	std::string failLog = "";
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}, {1.0, 1.0, 2.0}, \
													 {1.1, 1.1, 1.2},{0.0, 0.4, 1.0}, {1.5, 1.5, 0.0}	};
	std::vector<bettiBoundaryTableEntry> bettiTable;
	
	//Attempt to write stats
	// RET: bool
	if(!ws->writeStats("test","test")){ failLog += "writeOutput writeStats failed\n"; }
	
	//Remove the file
	std::remove("test_stats.csv");
	
	//Attempt to write CSV
	// RET: bool
	if(!ws->writeCSV("test","test")){ failLog += "writeOutput writeCSV(str,str) failed\n"; }
	
	//Attempt to write CSV
	// RET: bool
	if(!ws->writeCSV("test","test","test")){ failLog += "writeOutput writeCSV(str,str,str) failed\n"; }
	
	//Attempt to write CSV
	// RET: bool
	if(!ws->writeCSV(testValueArray,"test")){ failLog += "writeOutput writeCSV(vector<vector<double>>,str) failed\n"; }
	
	//Attempt to write CSV
	// RET: bool
	if(!ws->writeCSV(testValueArray,"test","test")){ failLog += "writeOutput writeCSV(vector<vector<double>>,str,str) failed\n"; }
	
	//Attempt to write MAT
	// RET: bool
	if(!ws->writeMAT(testValueArray,"test")){ failLog += "writeOutput writeMAT failed\n"; }
	
	//Attempt to write Barcodes
	// RET: bool
	if(!ws->writeBarcodes(bettiTable,"barcodeOutput")){ failLog += "writeOutput writeBarcodes(str,str) failed\n"; }
	
	//Remove the files
	std::remove("test.csv");
	std::remove("test.mat");
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: Write Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: Write Test Functions---------------------------\n";
	}
	
	return;
}

// TEST write Functions
void t_write_csv_functions(std::string &log){
	auto *ws = new writeOutput();
	auto *rs = new readInput();
	std::string failLog = "";
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}, {1.0, 1.0, 2.0}, \
													 {1.1, 1.1, 1.2},{0.0, 0.4, 1.0}, {1.5, 1.5, 0.0}	};
	
	//Attempt to write stats
	// RET: bool
	if(!ws->writeCSV(testValueArray,"testCSVOutput")){ failLog += "writeOutput writeCSV failed\n"; }
	
	auto testRetArray = rs->readCSV("testCSVOutput.csv");
	
	if(testValueArray != testRetArray)
		failLog += "writeOutput write_csv read error occurred";

	//Remove the file
	std::remove("testCSVOutput.csv");

	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: Write CSV Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: Write CSV Functions---------------------------\n";
	}

	return;
}


// TEST write Functions
void t_write_mat_functions(std::string &log){
	auto *ws = new writeOutput();
	auto *rs = new readInput();
	std::string failLog = "";
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}, {1.0, 1.0, 2.0}, \
													 {1.1, 1.1, 1.2},{0.0, 0.4, 1.0}, {1.5, 1.5, 0.0}	};
	
	//Attempt to write stats
	// RET: bool
	if(!ws->writeMAT(testValueArray,"testMATOutput")){ failLog += "writeOutput writeMAT failed\n"; }
	
	auto testRetArray = rs->readMAT("testMATOutput.mat");
	
	if(testValueArray != testRetArray)
		failLog += "writeOutput write_mat read error occurred";
	
	//Remove the file
	//std::remove("testMATOutput.mat");

	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: Write CSV Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: Write CSV Functions---------------------------\n";
	}

	return;
}

int main (int, char**){
	std::string log;
	t_write_functions(log);
	t_write_csv_functions(log);
	t_write_mat_functions(log);
	
	std::cout << std::endl << log << std::endl;
}
