#pragma once

// Header file for kMeansPlusPlus class - see kMeansPlusPlus.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class streamingUtils  {
  private:
	int num_clusters;			
	int num_iterations;			
  public:
	streamingUtils();
  void  kMeans(std::vector<std::vector<double>>& kHat);
  void ballKmeans(std::vector<std::vector<double>>& kHat);
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};

