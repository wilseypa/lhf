
/*
 * argParser hpp + cpp protoype and define a class for parsing command line
 * arguments and returning a map (dictionary) of parsed arguments
 *
 */

/**
 * @file argParser.hpp
 * @brief Defines the argParser class for parsing command-line arguments.
 */

#include <string>
#include <iostream>
#include "argParser.hpp"

std::map<std::string, std::string> argMap = {{"betaMesh","bmesh"},{"betaMode","bm"},{"beta","b"},{"alphaFilterationValue","afv"},{"nodeType","n"},{"reductionPercentage","rp"},{"maxSize","ms"},{"threads","t"},{"threshold","th"},{"scalar","s"},{"mpi","a"},{"mode","m"},{"dimensions","d"},{"iterations","r"},{"pipeline","p"},{"inputFile","i"},{"outputFile","o"},{"epsilon","e"},{"lambda","l"},{"debug","x"},{"complexType","c"},{"clusters","k"},{"preprocessor","pre"},{"upscale","u"},{"twist","w"},{"collapse","z"},{"seed","q"},{"involutedUpscale","iu"}, {"involuted","inv"}};
std::map<std::string, std::string> defaultMap = {{"betaMesh","null.csv"},{"betaMode","noMode"},{"beta","1"},{"alphaFilterationValue","50000"},{"nodeType","simplexNode"},{"reductionPercentage","10"},{"maxSize","2000"},{"threads","30"},{"threshold","250"},{"scalar","0.5"},{"mpi", "0"},{"mode", "standard"},{"dimensions","1"},{"iterations","250"},{"pipeline",""},{"inputFile","None"},{"outputFile","output"},{"epsilon","5"},{"lambda",".25"},{"debug","0"},{"complexType","simplexArrayList"},{"clusters","20"},{"preprocessor",""},{"upscale","false"},{"seed","-1"},{"twist","false"},{"collapse","false"},{"involutedUpscale","false"},{"involuted","false"}};
// argParse constructor, currently no needed information for the class constructor
/**
 * @brief argParser constructor.
 *
 * Currently, no additional information is needed for the class constructor.
 */
argParser::argParser(){

}

// argParser::defaultArguments ->

/**
 * @brief Sets the default values for command-line arguments.
 *
 * This function loops through each pair in argMap (maps shorthand arguments),
 * replacing the short argument if it exists, or mapping the default argument
 * if the long argument doesn't exist.
 *
 * @param map The map of argument names and their values to be modified.
 */

void argParser::defaultArguments(std::map<std::string, std::string>  &map){

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
}

/**
 * @brief Prints the usage information for the program.
 *
 * This function displays the command-line usage and available options
 * for the program, including their default values if applicable.
 */

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
	std::cout << "\t -l,--lambda <double>" << std::endl;
	std::cout << "\t\tDecay factor lambda for DenStream" << std::endl;
	std::cout << "\t\tdefault: .25" << std::endl;
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

/**
 * @brief Parses the command-line arguments and returns a map of arguments.
 *
 * This function takes the argument count and array of arguments from the
 * main function, processes the arguments, and returns a map (dictionary)
 * containing the parsed arguments along with their default values.
 *
 * @param argc Argument count from the main function.
 * @param argv Array of arguments from the main function.
 * @return A map containing the parsed arguments and their values.
 */

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

/**
 * @brief Prints the provided arguments to the console.
 *
 * This function takes a map (dictionary) of arguments and their values,
 * and prints them to the console in a user-friendly format.
 *
 * @param args A map containing the parsed arguments and their values.
 */

/**
 * @brief Prints the provided arguments to the console.
 *
 * This function takes a map (dictionary) of arguments and their values,
 * and prints them to the console in a user-friendly format.
 *
 * @param args A map containing the parsed arguments and their values.
 */

void argParser::printArguments(std::map<std::string,std::string> args){

	//Print the argument set
	std::cout << "Arguments being used: " << std::endl;
	for( const auto& sm_pair : args	){
		std::cout << "\t" << sm_pair.first << " (" << argMap[sm_pair.first] << ")   \t" << sm_pair.second << std::endl;
	}
	std::cout << std::endl;
	return;
}

/**
 * @brief Sets the pipeline for the given set of arguments.
 *
 * This function configures the pipeline for various modes and complex types
 * based on the provided arguments. It modifies the input map (dictionary) of
 * arguments to include the appropriate pipeline settings.
 *
 * @param args A map containing the parsed arguments and their values.
 */

void argParser::setPipeline(std::map<std::string, std::string>& args){
	//Configure the base pipeline from the complex storage type
	std::string basePipeline;


	// Handle iterative / involuted and upscaling flags
	if(args["mode"] == "iterUpscale" || args["involutedUpscale"] == "true"){
		args["upscale"] = "true";
	}

	
	// Handle complexType and corresponding simplexType	
	if(args["complexType"] == "simplexArrayList"){
		if(args["upscale"] == "true")
			basePipeline = "distMatrix.neighGraph.incrementalPersistence.upscale";
		else
			basePipeline = "distMatrix.neighGraph.incrementalPersistence";

	} else if (args["complexType"] == "simplexTree"){
		if(args["upscale"] == "true")
			basePipeline = "distMatrix.neighGraph.rips.fast.upscale";
		else
			basePipeline = "distMatrix.neighGraph.rips.fast";
			
	} else if(args["complexType"] == "alpha"){
		basePipeline = "distMatrix.alpha.fastPersistence";
		args["complexType"] = "alphaComplex";
		args["nodeType"] = "alphaNode";
		args["mode"] = "alpha";   
		args["pipeline"] = basePipeline;
		return;
	}else if(args["complexType"] == "graphInducedComplex" || args["complexType"] == "beta" || args["complexType"] == "betaGeneral" || args["complexType"] == "betaExtended"){
		args["nodeType"] = "alphaNode";
		args["mode"] = "alpha";
		if(args["complexType"] == "beta")
			basePipeline = "distMatrix.neighGraph.rips.fastPersistence";
		else if(args["complexType"] == "betaGeneral")
			basePipeline = "distMatrix.betaSkeletonBasedComplex.neighGraph.rips.fastPersistence";
		else if(args["complexType"] == "betaExtended")
			basePipeline = "distMatrix.betaSubSkeletonComplex.neighGraph.rips.fastPersistence";
		
		args["pipeline"] = basePipeline;
		args["complexType"] = "alphaComplex";
		return;
    } else if(args["complexType"] == "witness"){
		
		args["complexType"] = "witnessComplex";
		args["nodeType"] = "witnessNode";
		if(args["upscale"] == "true")
			basePipeline = "distMatrix.neighGraph.incrementalPersistence.upscale";
		else
			basePipeline = "distMatrix.neighGraph.incrementalPersistence";
			
    }
        
        
        
        
    if(args["mode"] == "mpi"){
		//Set up MPI pipeline ; requires setting pipeline, any complex storage
		if(args["pipeline"] == "")
			args["pipeline"] = basePipeline;

		//Requires preprocessor
		if(args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";

	}else if(args["mode"] == "standard"){
		//Set up standard mode;
		if(args["pipeline"] == "")
			args["pipeline"] = basePipeline;

	}else if(args["mode"] == "reduced"){
		//Set up reduced pipeline for centroid approximation

		//Requires preprocessor
		if(args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";

		if(args["pipeline"] == "")
			args["pipeline"] = basePipeline;

	} else if(args["mode"] == "upscale"){

		//Check if upscale was set or not, should be set in this mode
		if(args["upscale"] != "true"){
			args["upscale"] = "true";
			basePipeline += ".upscale";
		}

		//Requires preprocessor
		if(args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";
		if(args["pipeline"] == "")
			args["pipeline"] = basePipeline;

	} else if(args["mode"] == "stream"){
		if(args["preprocessor"] == "")
			args["preprocessor"] = "streamingkmeans";
		if(args["pipeline"] == "")
			args["pipeline"] = basePipeline;

	} else if(args["mode"] == "sw" || args["mode"] == "slidingwindow"){
		args["preprocessor"] = "";
		args["pipeline"] = "slidingwindow";
		args["complexType"] = "simplexTree";

	} else if(args["mode"] == "fast"){
		//Fast mode; uses the fast pipe regardless of simplexType

		if(args["upscale"] == "true") {
			args["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence.upscale";
		} else {
			args["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence";
		}

	} else if(args["mode"] == "naive" || args["mode"] == "naivewindow"){
		args["preprocessor"] = "";
		args["pipeline"] = "naivewindow";

	} else if(args["mode"] == "iterUpscale" || args["mode"] == "iter"){
		//Iterative upscaling pipe; requires preprocessor and uses basePipeline

		//Requires preprocessor
		if(args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";

		if(args["pipeline"] == "")
			args["pipeline"] = basePipeline;

	} else if(args["mode"] == "inc" || args["mode"] == "incremental"){
		//Incremental Mode; uses the incremental pipe regardless of simplexType

		if(args["complexType"] != "simplexArrayList")
			std::cout << "**SimplexArrayList is required for incremental persistence. Changing automatically.**" << std::endl;

		//Requires SAL
		args["complexType"] = "simplexArrayList";

		//Reset the pipe after complex change
		if(args["upscale"] == "true") {
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
		} else {
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence";
		}
	}

	return;
}
