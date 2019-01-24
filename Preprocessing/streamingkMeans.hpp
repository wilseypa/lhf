#pragma once

// Header file for streamingkMeans class - see streamingkMeans.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class streamingkMeans : public preprocessor {
  private:
	int num_clusters;			
	//int num_iterations;			
  public:
	streamingkMeans();
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};
