#pragma once

// Header file for kMeansPlusPlus class - see kMeansPlusPlus.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class kMeansPlusPlus : public preprocessor {
  private:
	int num_clusters;			
	int num_iterations;			
  public:
	kMeansPlusPlus();
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};

