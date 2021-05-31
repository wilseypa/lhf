#pragma once

// Header file for kMeansPlusPlus class - see kMeansPlusPlus.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

template<typename T>
class kMeansPlusPlus : public preprocessor<T> {
  private:
    int seed;
	int num_clusters;			
	int num_iterations;			
  public:
	kMeansPlusPlus();
    void runPreprocessor(pipePacket<T>&);
    bool configPreprocessor(std::map<std::string, std::string>&);
};

