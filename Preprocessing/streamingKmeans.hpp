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
  struct bucket { //bucket holds a coreset of "m" points
 //   int curSize;
    std::vector<std::vector<double>> points;
  //  struct point *points;
 //   struct point *spillover;  //describes points when a bucket reaches size m
  };


    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};
