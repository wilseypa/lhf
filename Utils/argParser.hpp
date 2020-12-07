#pragma once

// Header file for argParser class - see argParse.cpp for descriptions
#include <map>

class argParser {
  private:
  public:
    argParser();
    static void printUsage();
    static void printArguments(std::map<std::string,std::string>);
    static std::map<std::string, std::string> parse(int argc, char** argv);
    static void defaultArguments(std::map<std::string, std::string> &map);
    static void setPipeline(std::map<std::string, std::string>&);
};

