#pragma once

// Header file for streamingkMeans class - see streamingkMeans.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class streamingKmeans : public preprocessor {
  private:
	int num_clusters;			
	int num_iterations;		
  	
  public:
	streamingKmeans();
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};
