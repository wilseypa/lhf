#pragma once

// Header file for preprocessor class - see preprocessor.cpp for descriptions
#include <map>
#include "pipePacket.hpp"

class preprocessor {
  private:
  public:
	std::string procName;
	int debug;
    preprocessor();
    preprocessor* newPreprocessor(const std::string&);
    pipePacket runPreprocessorWrapper(pipePacket inData);
    virtual pipePacket runPreprocessor(pipePacket inData);
    virtual bool configPreprocessor(std::map<std::string, std::string> configMap);
};

