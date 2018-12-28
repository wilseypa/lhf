#pragma once

// Header file for argParser class - see argParse.cpp for descriptions
#include <map>

class argParser {
  private:
  public:
    argParser();
    void printUsage();
    std::map<std::string, std::string> parse(int argc, char** argv);
    std::map<std::string, std::string> defaultArguments(std::map<std::string, std::string>  &map);
};

