
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

struct Args
{
	std::string shorthand;	  // Shorthand name for the argument
	std::string defaultValue; // Default value for the argument
	std::string description;  // Explanation or description of the argument
	std::string d_type;		  // DataType of the Argument

	Args() : shorthand(""), defaultValue(""), description(""), d_type("") {}

	Args(const std::string &sh, const std::string &def, const std::string &desc, const std::string &d_t)
		: shorthand(sh), defaultValue(def), description(desc), d_type(d_t)
	{
		if (description == "")
			description = "No Description available";
	}
};

static std::map<std::string, Args> argMap = {
	{"betaMesh", Args("bmesh", "null.csv", "", "")},
	{"betaMode", Args("bm", "noMode", "", "")},
	{"beta", Args("b", "1", "", "<int>")},
	{"alphaFilterationValue", Args("afv", "50000", "", "")},
	{"nodeType", Args("n", "simplexNode", "", "")},
	{"reductionPercentage", Args("rp", "10", "", "")},
	{"maxSize", Args("ms", "2000", "", "")},
	{"threads", Args("t", "30", "", "<int>")},
	{"threshold", Args("th", "250", "", "")},
	{"scalar", Args("s", "0.5", "", "")},
	{"mpi", Args("a", "0", "", "<int>")},
	{"mode", Args("m", "standard", "Sets the mode for LHF to run in", "(standard|reduced|upscale|sw)")},
	{"dimensions", Args("d", "1", "Sets the maximum homology dimension to compute (H_d)", "<int>")},
	{"iterations", Args("r", "250", "", "<int>")},
	{"pipeline", Args("p", "", "", "")},
	{"inputFile", Args("i", "None", "Filename (csv) for LHF input", "<filename>")},
	{"outputFile", Args("o", "output", "Filename for LHF output", "<filename>")},
	{"epsilon", Args("e", "5", "Maximum epsilon threshold", "<float>")},
	{"lambda", Args("l", ".25", "Decay factor lambda for DenStream", "")},
	{"debug", Args("x", "0", "", "<int(0|1)>")},
	{"complexType", Args("c", "simplexArrayList", "", "")},
	{"clusters", Args("k", "20", "", "<int>")},
	{"preprocessor", Args("pre", "", "", "")},
	{"upscale", Args("u", "false", "", "<bool>")},
	{"seed", Args("q", "-1", "", "")},
	{"twist", Args("w", "false", "", "<bool>")},
	{"collapse", Args("z", "false", "", "<bool>")},
	{"involutedUpscale", Args("iu", "false", "", "<bool>")},
	{"involuted", Args("inv", "false", "", "<bool>")}};

/**
 * @brief Constructor for argParser class.
 *
 * Currently, no additional information is needed for the class constructor.
 */

argParser::argParser()
{
}

/**
 * @brief Sets the default values for command-line arguments.
 *
 * This function loops through each pair in argMap (maps shorthand arguments),
 * replacing the short argument if it exists, or mapping the default argument
 * if the long argument doesn't exist.
 *
 * @param map The map of argument names and their values to be modified.
 */

void argParser::defaultArguments(std::map<std::string, std::string> &map)
{

	// Loop through each pair in argMap (maps shorthand arguments)
	//	if the short argument exists, replace it
	//	else, if the long argument doesn't exist map the default argument
	for (const auto &arg_pair : argMap)
	{
		auto i = map.find(arg_pair.second.shorthand);
		if (i != map.end())
		{
			map[arg_pair.first] = i->second;
			map.erase(i);
		}
		else
		{
			i = map.find(arg_pair.first);
			if (i == map.end())
			{
				map[arg_pair.first] = arg_pair.second.defaultValue;
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

void argParser::printUsage()
{
	std::cout << "Usage: " << std::endl;
	std::cout << std::endl;
	std::cout << "\t ./LHF -i <inputfile> -o <outputfile>" << std::endl;
	std::cout << std::endl;
	std::cout << "Additional Options:" << std::endl;
	for (const auto &arg : argMap)
	{
		std::cout << "\t -" << arg.second.shorthand << ",--" << arg.first << " " << arg.second.d_type << std::endl;
		std::cout << "\t\t" << arg.second.description << std::endl;
		std::cout << "\t\t\tDefault: " << arg.second.defaultValue << std::endl;
		std::cout << std::endl;
	}

	return;
}

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

std::map<std::string, std::string> argParser::parse(int argc, char **argv)
{
	std::map<std::string, std::string> retVal;
	// Remove and map the initial arguments (may be led with -- or -)
	for (int i = 1; i < argc; i = i + 2)
	{
		if (std::string(argv[i]).substr(0, 2) == "--")
		{
			retVal[std::string(argv[i]).substr(2)] = argv[i + 1];
		}
		else if (std::string(argv[i]).substr(0, 1) == "-")
		{
			retVal[std::string(argv[i]).substr(1)] = argv[i + 1];
		}
		else
		{
			std::cout << "Invalid argument: " << argv[i] << std::endl
					  << std::endl;
			printUsage();
		}
	}

	// Call defaultArguments to translate shorthand arguments and add defaults to map
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

void argParser::printArguments(const std::map<std::string, std::string> &args)
{

	// Print the argument set
	std::cout << "Arguments being used: " << std::endl;
	for (const auto &sm_pair : args)
	{
		std::cout << "\t" << sm_pair.first << " (" << argMap[sm_pair.first].shorthand << ")   \t" << sm_pair.second << std::endl;
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

void argParser::setPipeline(std::map<std::string, std::string> &args)
{

	// Handle iterative / involuted and upscaling flags
	if (args["mode"] == "iterUpscale" || args["involutedUpscale"] == "true")
	{
		args["upscale"] = "true";
	}

	// Handle witness complexes
	if (args["complexType"] == "witness" || args["complexType"] == "witnessComplex")
	{
		args["complexType"] = "witnessComplex";
		args["nodeType"] = "witnessNode";
	}

	/**
	 * MPI MODE:
	 *  VR with preprocessor, partitions farmed to MPI nodes
	 */
	if (args["mode"] == "mpi")
	{
		args["mpi"] = 1;

		// Set up MPI pipeline ; requires setting pipeline, any complex storage
		if (args["pipeline"] == "")
		{
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence";

			if (args["upscale"] == "true")
				args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
		}

		// Requires preprocessor
		if (args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";
	}

	/**
	 * STANDARD MODE:
	 * 	VR with incremental
	 */
	else if (args["mode"] == "standard")
	{
		// Set up standard mode;
		if (args["pipeline"] == "")
		{
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence";

			if (args["upscale"] == "true")
				args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
		}
	}

	/**
	 * ALPHA MODE:
	 * 	Triangulation with Alpha (gabriel) filtration
	 */
	else if (args["mode"] == "alpha" || args["complexType"] == "alpha" || args["complexType"] == "alphaComplex")
	{

		// TODO: Nick come check these are still correct after splitting betaComplex
		args["nodeType"] = "alphaNode";
		args["mode"] = "alpha";
		args["complexType"] = "alphaComplex";
		args["pipeline"] = "distMatrix.alpha.fastPersistence";

		if (args["upscale"] == "true")
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
	}

	/**
	 * BETA MODE:
	 * 	Triangulation with Beta-sparsification
	 */
	else if (args["mode"] == "beta" || args["complexType"] == "graphInducedComplex" || args["complexType"] == "betaGeneral" || args["complexType"] == "betaExtended")
	{
		args["nodeType"] = "alphaNode";
		args["mode"] = "beta";

		if (args["complexType"] == "betaGeneral")
			args["pipeline"] = "distMatrix.betaSkeletonBasedComplex.neighGraph.rips.fastPersistence";
		else if (args["complexType"] == "betaExtended")
			args["pipeline"] = "distMatrix.betaSubSkeletonComplex.neighGraph.rips.fastPersistence";

		args["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence";
		args["complexType"] = "betaComplex";
	}

	/**
	 * WEIGHTED DELAUNAY MODE:
	 * 	Triangulation with Alpha (gabriel) filtration
	 */
	else if (args["mode"] == "wd" || args["mode"] == "weightedAlpha")
	{

		// TODO: Nick come check these are still correct after splitting betaComplex
		args["nodeType"] = "alphaNode";
		args["mode"] = "weightedAlpha";
		args["complexType"] = "alphaComplex";
		args["pipeline"] = "distMatrix.qhull.incrementalPersistence";

		if (args["upscale"] == "true")
			args["pipeline"] = "distMatrix.qhull.incrementalPersistence.upscale";
	}

	/**
	 * REDUCED MODE:
	 * 	VR with preprocessor
	 */
	else if (args["mode"] == "reduced")
	{
		// Set up reduced pipeline for centroid approximation

		// Requires preprocessor
		if (args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";

		if (args["pipeline"] == "")
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence";

		/**
		 * UPSCALE MODE:
		 * 	VR with preprocessor, upscaling
		 */
	}
	else if (args["mode"] == "upscale")
	{

		args["upscale"] = "true";

		// Requires preprocessor
		if (args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";

		if (args["pipeline"] == "")
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
	}

	/**
	 * STREAM MODE:
	 * 	Streaming VR with preprocessor
	 */
	else if (args["mode"] == "stream")
	{
		if (args["preprocessor"] == "")
			args["preprocessor"] = "streamingkmeans";
		if (args["pipeline"] == "")
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence";
		;
	}

	/**
	 * SLIDING WINDOW MODE:
	 * 	VR with decision-based sliding window
	 */
	else if (args["mode"] == "sw" || args["mode"] == "slidingwindow")
	{
		args["preprocessor"] = "";
		args["pipeline"] = "slidingwindow";
		args["complexType"] = "simplexTree";
	}

	/**
	 * FAST MODE:
	 * 	VR without fast cofacet indexing
	 */
	else if (args["mode"] == "fast")
	{
		if (args["upscale"] == "true")
		{
			args["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence.upscale";
		}
		else
		{
			args["pipeline"] = "distMatrix.neighGraph.rips.fastPersistence";
		}
	}

	/**
	 * NAIVE WINDOW MODE:
	 * 	Simple sliding window VR
	 */
	else if (args["mode"] == "naive" || args["mode"] == "naivewindow")
	{
		args["preprocessor"] = "";
		args["pipeline"] = "naivewindow";
	}

	/**
	 * ITERATIVE UPSCALING MODE:
	 * 	VR with preprocessor, iterative upscaling
	 */
	else if (args["mode"] == "iterUpscale" || args["mode"] == "iter")
	{
		// Iterative upscaling pipe; requires preprocessor and uses basePipeline

		// Requires preprocessor
		if (args["preprocessor"] == "")
			args["preprocessor"] = "kmeans++";

		if (args["pipeline"] == "")
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
	}

	/**
	 * INCREMENTAL MODE:
	 * 	VR, same as standard mode
	 */
	else if (args["mode"] == "inc" || args["mode"] == "incremental")
	{
		// Incremental Mode; uses the incremental pipe regardless of simplexType

		if (args["complexType"] != "simplexArrayList")
			std::cout << "**SimplexArrayList is required for incremental persistence. Changing automatically.**" << std::endl;

		// Requires SAL
		args["complexType"] = "simplexArrayList";

		// Reset the pipe after complex change
		if (args["upscale"] == "true")
		{
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence.upscale";
		}
		else
		{
			args["pipeline"] = "distMatrix.neighGraph.incrementalPersistence";
		}
	}

	return;
}
