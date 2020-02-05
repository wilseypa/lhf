#pragma once

// Header file for densityUtils class - see densityUtils.cpp for descriptions
#include <map>
#include <vector>
#include "preprocessor.hpp"
#include "kdTree.hpp"

class dbscan : public preprocessor {
  private:
    int minPoints;
    double epsilon;
  	
    void expandCluster(const std::vector<std::vector<double>> &data, 
                       std::vector<int> &labels, 
                       std::vector<size_t> &neighbors, 
                       int clusterLabel,
                       kdTree &tree);

  public:
    dbscan();
    std::vector<int> cluster(const std::vector<std::vector<double>> &data);
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
}; 
