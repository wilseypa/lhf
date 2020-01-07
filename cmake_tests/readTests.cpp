#include <iostream>
#include <vector>
#include "readInput.hpp"

// TEST read Functions
void t_read_functions(std::string &log){
	auto *rs = new readInput();
	std::string failLog = "";
	std::vector<double> a = {0};
	
	//Attempt to read csv
	// RET: vector<vector<double>>
	if(rs->readCSV("test").size() != 0){ failLog += "readInput readCSV failed\n"; }
	
	//Attempt to read mat
	// RET: vector<vector<double>>
	if(rs->readMAT("test").size() != 0){ failLog += "readInput readCSV failed\n"; }
	
	//Attempt to initialize stream
	// RET: bool
	if(rs->streamInit("test")){ failLog += "readInput streamInit failed\n"; }
	
	//Attempt to read from stream
	// RET: bool
	if(rs->streamRead(a)){ failLog += "readInput streamInit failed\n"; }
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: Read Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: Read Test Functions---------------------------\n";
	}
	
	return;
}

int main (int, char**){
	std::string log;
	t_read_functions(log);
	
	std::cout << std::endl << log << std::endl;
}
