#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include "readInput.hpp"
#include "argParser.hpp"
#include "basePipe.hpp"
#include "writeOutput.hpp"
#include "pipePacket.hpp"

using namespace std;

int main(int argc, char* argv[]){
    auto *rs = new readInput();
    auto *ap = new argParser();
    auto *ws = new writeOutput();
    auto *wD = new pipePacket();
    
    auto args = ap->parse(argc, argv);
	
	
	//Read data from inputFile CSV
    wD->workData.originalData = rs->readCSV(args["inputFile"]);


	if(wD->workData.originalData.size() > 0){
		
		//Add data to our pipePacket
		wD->workData.workingData = wD->workData.originalData;
	
		// Begin processing parts of the pipeline
		// DataInput -> A -> B -> ... -> DataOutput
		// Parsed by "." -> i.e. A.B.C.D
		auto pipe = args.find("pipeline");
		if(pipe != args.end()){
			auto pipeFuncts = std::string(args["pipeline"]);
			auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;
			
			//For each '.' separated pipeline function (count of '.' + 1 -> lim)
			for(unsigned i = 0; i < lim; i++){
				auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
				pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);
				
				//Build the pipe component, configure and run
				auto *bp = new basePipe();
				auto *cp = bp->newPipe(curFunct);
				
				if(cp->configPipe(args)){
					*wD = cp->runPipe(*wD);
				} else {
					cout << "Failed to configure pipeline: " << args["pipeline"] << endl;
				}
			}
		}
		else {
			cout << "Failed to find a suitable pipeline, exiting..." << endl;
			return 1;
		}
		
		pipe = args.find("outputFile");
		/*if(pipe != args.end()){
			if (args["outputFile"] == "console"){
				ws->writeConsole(wD);
			
			}
		}*/
		
	}
    return 0;
}
