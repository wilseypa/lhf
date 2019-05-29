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
double streamingKmeans:: dotProd(const std::vector<double>& a, const std::vector<double>& b);

double streamingKmeans::dotProd2D(std::vector<std::vector<double>>& a, std::vector<std::vector<double>> & b);

void streamingKmeans:: approxNearestNeighbor(std::vector<std::vector<double>> approxFacilities, std::vector<double> omega, int x, int size, pipePacket(inData)));

int streamingKmeans:: binarySearch(std::vector<std::vector<double>> approxFacilities, std::vector<double> omega, int n, double target);

double streamingKmeans::randDouble(); // Returns a random double in the range [0,1)
bool streamingKmeans::prob(double f); 
int streamingKmeans::random(int low, int high);   // returns random int in the range [low, high)

    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
}; 
