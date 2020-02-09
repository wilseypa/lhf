 
/*
 * argParser hpp + cpp protoype and define a class for parsing command line
 * arguments and returning a map (dictionary) of parsed arguments
 * 
 */

#include <string>
#include <iostream>
#include "argParser.hpp"


std::map<std::string, std::string> argMap = { {"mode","m"},{"dimensions","d"},{"iterations","r"},{"pipeline","p"},{"inputFile","i"},{"outputFile","o"},{"epsilon","e"},{"debug","x"},{"complexType","c"},{"clusters","k"},{"preprocessor","pre"},{"upscale","u"},{"twist","t"},{"collapse","z"}};
std::map<std::string, std::string> defaultMap = { {"mode", "standard"},{"dimensions","1"},{"iterations","250"},{"pipeline",""},{"inputFile","None"},{"outputFile","output.csv"},{"epsilon","5"},{"debug","0"},{"complexType","simplexTree"},{"clusters","5"},{"preprocessor",""},{"upscale","false"},{"twist","false"},{"collapse","false"}};

// argParse constructor, currently no needed information for the class constructor
argParser::argParser(){

}

// argParser::defaultArguments -> 
std::map<std::string, std::string> argParser::defaultArguments(std::map<std::string, std::string>  &map){
	
	//Loop through each pair in argMap (maps shorthand arguments)
	//	if the short argument exists, replace it
	//	else, if the long argument doesn't exist map the default argument
	for(const auto& arg_pair : argMap) {
		auto i = map.find(arg_pair.second);
		if(i != map.end()){
			map[arg_pair.first] = map[arg_pair.second];
			map.erase(i);
		} else {
			i = map.find(arg_pair.first);
			if(i == map.end()){
				map[arg_pair.first] = defaultMap[arg_pair.first];
			}
		}
	}
	
	return map;
}


void argParser::printUsage(){
	std::cout << "Usage: " << std::endl;
	std::cout << std::endl;
	std::cout << "\t ./LHF -i <inputfile> -o <outputfile>" << std::endl;
	std::cout << std::endl;
	std::cout << "Additional Options:" << std::endl;
	std::cout << "\t -i,--inputFile <filename>" << std::endl;
	std::cout << "\t\tFilename (csv) for LHF input"  << std::endl;
	std::cout << std::endl; 
	std::cout << "\t -o,--outputFile <filename>" << std::endl;
	std::cout << "\t\tFilename for LHF output" << std::endl;
	std::cout << std::endl; 
	std::cout << "\t -e,--epsilon <int>" << std::endl;
	std::cout << "\t\tMaximum epsilon threshold" << std::endl;
	std::cout << "\t\tdefault: 5" << std::endl;
	std::cout << std::endl;
	std::cout << "\t -m,--mode (standard|reduced|upscale|sw)" << std::endl;
	std::cout << "\t\tSets the mode for LHF to run in" << std::endl;
	std::cout << "\t\t\tdefault: standard" << std::endl;
	std::cout << std::endl;
	std::cout << "\t -d,--dimensions <int>" << std::endl;
	std::cout << "\t\tSets the maximum homology dimension to compute (H_d)" << std::endl;
	std::cout << "\t\t\tdefault: 1" << std::endl;
	std::cout << std::endl;
	
	return;
}


// argParser::parse -> Parse through arguments and return a map (dictionary)
//		-argc - argument count from Main()
//		-argv - array of arguments from Main()
// 
std::map<std::string, std::string> argParser::parse(int argc, char** argv){
	std::map<std::string,std::string> retVal;
	//Remove and map the initial arguments (may be led with -- or -)
	for(int i = 1; i < argc; i=i+2){
		if (std::string(argv[i]).substr(0,2) == "--"){
			retVal[std::string(argv[i]).substr(2)] = argv[i+1];
		} else if (std::string(argv[i]).substr(0,1) == "-"){
			retVal[std::string(argv[i]).substr(1)] = argv[i+1];
		} else {
			std::cout << "Invalid argument: " << argv[i] << std::endl << std::endl;
			printUsage();
		}
	}
	
	//Call defaultArguments to translate shorthand arguments and add defaults to map
	defaultArguments(retVal);

	return retVal;
}

void argParser::printArguments(std::map<std::string,std::string> args){
	
	//Print the argument set
	std::cout << "Arguments being used: " << std::endl;
	for( const auto& sm_pair : args	){
		std::cout << "\t" << sm_pair.first << " (" << argMap[sm_pair.first] << ")   \t" << sm_pair.second << std::endl;
	}
	std::cout << std::endl;
	return;
}

void argParser::setPipeline(std::map<std::string, std::string>& args){
	if(args["pipeline"] == ""){
		if(args["mode"] == "standard")
			args["pipeline"] = "distMatrix.neighGraph.rips.persistence";
		else if(args["mode"] == "reduced"){
			if(args["preprocessor"] == "")
				args["preprocessor"] = "kmeans++";
			args["pipeline"] = "distMatrix.neighGraph.rips.persistence";
		} else if(args["mode"] == "upscale"){
			if(args["preprocessor"] == "")
				args["preprocessor"] = "kmeans++";
			args["upscale"] = "true";
			args["pipeline"] = "distMatrix.neighGraph.rips.persistence.boundary.upscale";
		} else if(args["mode"] == "stream"){
			if(args["preprocessor"] == "")
				args["preprocessor"] = "streamingkmeans";
			args["pipeline"] = "distMatrix.neighGraph.rips.persistence";
		} else if(args["mode"] == "sw" || args["mode"] == "slidingwindow"){
			args["preprocessor"] = "";
			args["pipeline"] = "slidingwindow";
			args["upscale"] = "false";
			args["complexType"] = "simplexTree";
		} else if(args["mode"] == "fast"){
			args["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence";
			args["upscale"] = "false";
			args["complexType"] = "simplexTree";
		}	
	}
	
	return;
}


