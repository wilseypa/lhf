 
/*
 * argParser hpp + cpp protoype and define a class for parsing command line
 * arguments and returning a map (dictionary) of parsed arguments
 * 
 */

#include <string>
#include <iostream>
#include "argParser.hpp"


std::map<std::string, std::string> argMap = { {"dimensions","d"},{"iterations","r"},{"pipeline","p"},{"inputFile","i"},{"outputFile","o"},{"epsilon","e"} };
std::map<std::string, std::string> defaultMap = { {"dimensions","3"},{"iterations","1000"},{"pipeline","default"},{"inputFile","None"},{"outputFile","console"},{"epsilon","0.5"} };

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
				std::cout << "Invalid argument: " << argv[i] << std::endl;
		}
	}
	
	//Call defaultArguments to translate shorthand arguments and add defaults to map
	argParser::defaultArguments(retVal);
	
	//Print the argument set
	for( const auto& sm_pair : retVal){
		std::cout << sm_pair.first << '\t' << sm_pair.second << std::endl;
	}
	std::cout << std::endl;
	return retVal;
}




