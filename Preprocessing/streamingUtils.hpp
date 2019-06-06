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
  std::vector<std::vector<double>> kMeans(std::vector<std::vector<double>>& kHat);
  std::vector<std::vector<double>> ballKmeans(std::vector<std::vector<double>>& kHat);

  pipePacket runPreprocessor(pipePacket inData);
  bool configPreprocessor(std::map<std::string, std::string> configMap);
};

