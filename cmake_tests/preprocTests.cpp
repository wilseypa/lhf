#include <iostream>
#include <vector>
#include "preprocessor.hpp"
#include "pipePacket.hpp"

// TEST preprocessor Functions
void t_preproc_functions(std::string &log){
	/*
	std::string failLog = "";
	auto *pack = new pipePacket("simplexArrayList", 5.0, 2);
	preprocessor *testPreproc = new preprocessor();
	
	std::vector<unsigned> testValue;
	std::vector<std::vector<double>> testValueArray;
	std::map<std::string, std::string> testConfig;
	
	//Attempt to output preprocessor data
	// RET: void 
	testPreproc->outputData(testValue);
	
	//Attempt to output preprocessor data
	// RET: void
    testPreproc->outputData(testValueArray);
    
    //Configure the pipe with arguments
    // RET: bool
    testPreproc->configPreprocessor(testConfig);
	
	//Run the preprocessor selected preprocessor
    testPreproc->runPreprocessor(*pack);
    
	//Run the pipe wrapper selected preprocessor
    testPreproc->runPreprocessorWrapper(*pack);
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: preprocessor Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: preprocessor Test Functions---------------------------\n";
	}
	
	* */
	return;
}


// TEST simplexArrayList Functions
void t_preproc_empty_functions(std::string &log, std::string type){
	/*
	std::string failLog = "";
	auto *pack = new pipePacket("simplexArrayList", 5.0, 2);
	preprocessor *testPreproc = preprocessor::newPreprocessor(type);
	std::vector<unsigned> testValue = {0, 1, 2};
	std::vector<std::vector<double>> testValueArray {{0.0, 1.0, 2.0},{2.0, 1.0, 0.0}};
	
	std::map<std::string, std::string> testConfig = {{"epsilon","5.0"}};
	
	//Attempt to output preprocessor data
	// RET: void 
	testPreproc->outputData(testValue);
	
	//Attempt to output preprocessor data
	// RET: void
    testPreproc->outputData(testValueArray);
    
    //Configure the pipe with arguments
    // RET: bool
    testPreproc->configPreprocessor(testConfig);
	
	//Run the preprocessor selected preprocessor
    testPreproc->runPreprocessor(*pack);
    
	//Run the pipe wrapper selected preprocessor
    testPreproc->runPreprocessorWrapper(*pack);
    
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: " + type + " Preprocessor Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: " + type + " Preprocessor Test Functions---------------------------\n";
	}
	
	* */
	return;
}

int main (int, char**){
	std::string log;
	t_preproc_functions(log);
	
	
	for(std::string type : {"kmeans++","streamingKmeans"}){ 
		try{t_preproc_empty_functions(log, type);}
		catch(const std::exception){log += "FAILED: " + type + " Empty Test Functions---------------------------\n";}			
	}
	
	std::cout << std::endl << std::endl << log << std::endl;
}
