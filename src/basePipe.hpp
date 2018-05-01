#pragma once

// Header file for basePipe class - see basePipe.cpp for descriptions
#include <map>

class basePipe {
  private:
  public:
	std::string pipeType;
    basePipe();
    basePipe* newPipe(const std::string&);
    virtual std::vector<std::vector<double>> runPipe(std::vector<std::vector<double>> inData);
    virtual bool configPipe(std::map<std::string, std::string> configMap);
};

