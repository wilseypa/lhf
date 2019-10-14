#pragma once

// Header file for densityUtils class - see densityUtils.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class dbscan : public preprocessor {
  private:
	
  	
  public:
  dbscan();
//  std::vector<std::vector<double>> cluster
      std::vector<int> cluster(std::vector<std::vector<double>> &data);
  int calcClusterVec(std::vector<std::vector<double>>& data, int i, double epsilon);
  int calcClusterPipePacket(pipePacket inData, int i, double epsilon);
  int expandCluster(std::vector<std::vector<double>>& data, std::vector<std::vector<double>>& neighbors, std::vector<int>& neighborTracker, std::vector<int>& clusterLabel, int clusterIndex, double epsilon, int minPoints, std::vector<int>& processed);

  std::vector<std::vector<double>> neighborQuery(std::vector<std::vector<double>>& data, std::vector<int>& neighborLabel, int i, double epsilon);
  pipePacket runPreprocessor(pipePacket inData);
  bool configPreprocessor(std::map<std::string, std::string> configMap);
}; 
