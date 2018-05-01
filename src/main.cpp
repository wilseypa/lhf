#include <iostream>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include "readInput.hpp"
#include "argParser.hpp"
#include "basePipe.hpp"

using namespace std;

int main(int argc, char* argv[]){
    auto *rs = new readInput();
    auto *ap = new argParser();
    
    auto args = ap->parse(argc, argv);
	
	
	
    auto result = rs->readCSV(args["inputFile"]);


	if(result.size() > 0){
		cout << "_________CSV RESULTS__________" << endl;
		
		for(unsigned i = 0; i < result.size(); i++){
			for(unsigned j = 0; j < result[i].size(); j++){
				cout << result[i][j] << "\t";
			}
			cout << endl;
		}
	
		// Begin processing parts of the pipeline
		// DataInput -> A -> B -> ... -> DataOutput
		// Parsed by "." -> i.e. DataInput.A.B.DataOutput
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
				bp = bp->newPipe(curFunct);
				cout << typeid(bp).name() << endl;
				if(bp->configPipe(args)){
					result = bp->runPipe(result);
				} else {
					cout << "Failed to configure pipeline: " << args["pipeline"] << endl;
				}
			}
		}
		else {
			cout << "Failed to find a suitable pipeline, exiting..." << endl;
			return 1;
		}
	}
    return 0;
}
