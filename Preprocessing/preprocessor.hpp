#pragma once

// Header file for preprocessor class - see preprocessor.cpp for descriptions
#include <map>
#include "pipePacket.hpp"

class preprocessor {
  private:
  public:
	bool configured = false;
	std::string procName = "preprocessor";
	bool debug = 0;
	std::string outputFile;
	utils ut;
    preprocessor();
    virtual ~preprocessor(){};
    static preprocessor* newPreprocessor(const std::string&);
    void runPreprocessorWrapper(pipePacket &inData);
    virtual void runPreprocessor(pipePacket &inData);
    virtual void outputData(pipePacket&);
	virtual void outputData(std::vector<unsigned>);
	virtual void outputData(std::vector<std::vector<double>>);
    virtual bool configPreprocessor(std::map<std::string, std::string> &configMap);
};

