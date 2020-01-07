#include <iostream>
#include <vector>
#include "writeOutput.hpp"

// TEST write Functions
void t_write_functions(std::string &log){
	auto *ws = new writeOutput();
	std::string failLog = "";
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}, {1.0, 1.0, 2.0}, \
													 {1.1, 1.1, 1.2},{0.0, 0.4, 1.0}, {1.5, 1.5, 0.0}	};
	
	//Attempt to write stats
	// RET: bool
	if(!ws->writeStats("test","test")){ failLog += "writeOutput writeStats failed\n"; }
	
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
	if(!ws->writeBarcodes(testValueArray,"test")){ failLog += "writeOutput writeBarcodes failed\n"; }
	
	//Attempt to write Barcodes
	// RET: bool
	if(!ws->writeBarcodes("test","test")){ failLog += "writeOutput writeBarcodes(str,str) failed\n"; }
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: Write Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: Write Test Functions---------------------------\n";
	}
	std::cout << "Exiting write function testing..." << std::endl;
	
	return;
}

int main (int, char**){
	std::string log;
	t_write_functions(log);
	
	std::cout << std::endl << log << std::endl;
}
