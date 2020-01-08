#pragma once

// Header file for preprocessor class - see preprocessor.cpp for descriptions
#include <map>
#include "pipePacket.hpp"

class preprocessor {
  private:
  public:
	bool configured = false;
	std::string procName = "preprocessor";
	int debug;
	std::string outputFile;
	utils ut;
    preprocessor();
    preprocessor* newPreprocessor(const std::string&);
    pipePacket runPreprocessorWrapper(pipePacket inData);
    virtual pipePacket runPreprocessor(pipePacket inData);
	virtual void outputData(std::vector<unsigned>);
	virtual void outputData(std::vector<std::vector<double>>);
    virtual bool configPreprocessor(std::map<std::string, std::string> configMap);
};

