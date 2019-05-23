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

void streamingKmeans:: approxNearestNeighbor(std::vector<std::vector<double>> facilities, float dotProd, int n, float distSquare, int dim );
void streamingKmeans:: binarySearch(std::vector<std::vector<double>> approxFacilities, int n, double target);
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
}; 
