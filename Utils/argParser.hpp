#pragma once

// Header file for argParser class - see argParse.cpp for descriptions
#include <map>

class argParser {
  private:
  public:
    argParser();
    void printUsage();
	void printArguments(std::map<std::string,std::string>);
	
	/** argParser::parse(int, char**)
	 * 
	 *  @brief Parse through arguments and return a map (dictionary)
	 *	
	 *  @param argc argument count from Main()
	 *	@param argv array of arguments from Main()
	*/ 
    std::map<std::string, std::string> parse(int argc, char** argv);
    std::map<std::string, std::string> defaultArguments(std::map<std::string, std::string>  &map);
    void setPipeline(std::map<std::string, std::string>&);
};

