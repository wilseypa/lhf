#pragma once

// Header file for kMeansPlusPlus class - see kMeansPlusPlus.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class kMeansPlusPlus : public preprocessor {
  private:
    int seed;
	int num_clusters;			
	int num_iterations;			
  public:
	kMeansPlusPlus();
    void runPreprocessor(pipePacket&);
    bool configPreprocessor(std::map<std::string, std::string>&);
};

