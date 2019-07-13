#pragma once

// Header file for densityUtils class - see densityUtils.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class densityUtils : public preprocessor {
  private:
		
  	
  public:
  densityUtils();
  std::vector<std::vector<double>> dbscan(std::vector<std::vector<double>>& data);
  pipePacket runPreprocessor(pipePacket inData);
  bool configPreprocessor(std::map<std::string, std::string> configMap);
}; 
