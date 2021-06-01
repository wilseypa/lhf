#pragma once

// Header file for kMeansPlusPlus class - see kMeansPlusPlus.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

template<typename nodeType>
class kMeansPlusPlus : public preprocessor<nodeType> {
  private:
    int seed;
	int num_clusters;			
	int num_iterations;			
  public:
	kMeansPlusPlus();
    void runPreprocessor(pipePacket<nodeType>&);
    bool configPreprocessor(std::map<std::string, std::string>&);
};

