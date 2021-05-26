#include <iostream>
#include <vector>
#include "basePipe.hpp"
#include "pipePacket.hpp"

// TEST basePipe Functions
void t_pipe_functions(std::string &log){
	/*std::string failLog = "";
	auto *pack = new pipePacket("simplexArrayList", 5.0, 2);
	basePipe *testPipe = new basePipe();
	
	std::map<std::string, std::string> testConfig;
	
	//Attempt to output pipe data without initializing pipe
	// RET: void 
	testPipe->outputData(*pack);
	
	//Attempt to run the pipe without initializing
	// RET: void
    testPipe->runPipe(*pack);
    
    //Configure the pipe with arguments
    // RET: bool
    if(!testPipe->configPipe(testConfig)){ failLog += "basePipe config failed\n"; }
	
	//Run the pipe wrapper with no selected pipe
    testPipe->runPipeWrapper(*pack);
	
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: basePipe Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: basePipe Test Functions---------------------------\n";
	}*/
	return;
}


// TEST simplexArrayList Functions
void t_pipe_empty_functions(std::string &log, std::string type, std::string complexType){
	/*std::string failLog = "";
	auto *pack = new pipePacket(complexType, 5.0, 2);
	basePipe *testPipe = basePipe::newPipe(type, complexType);
	
	std::map<std::string, std::string> testConfig = {{"epsilon","5.0"}};
	
	//Attempt to output pipe data without initializing pipe
	// RET: void 
	testPipe->outputData(*pack);
	
	//Attempt to run the pipe without initializing
	// RET: void
    testPipe->runPipe(*pack);
    
    //Configure the pipe with arguments
    // RET: bool
    testPipe->configPipe(testConfig);
	
	//Run the pipe wrapper with no selected pipe
    testPipe->runPipeWrapper(*pack);
	//Output log status to calling function
	if(failLog.size() > 0){
		log += "FAILED: " + type + " Empty Test Functions---------------------------\n" + failLog;	
	} else {
		 log += "PASSED: " + type + " Empty Test Functions---------------------------\n";
	}
	*/
	return;
}

int main (int, char**){
	std::string log;
	t_pipe_functions(log);
	
	
	for(std::string type : {"distMatrix","neighGraph","upscale","persistence","slidingwindow","fastPersistence"}){ 
		try{t_pipe_empty_functions(log, type, "simplexArrayList");}
		catch(const std::exception){log += "FAILED: " + type + " Empty Test Functions---------------------------\n";}			
	}
	
	std::cout << std::endl << std::endl << log << std::endl;
}
