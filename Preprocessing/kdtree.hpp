#pragma once

// Header file for kMeansPlusPlus class - see kMeansPlusPlus.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class kdtree : public preprocessor
{
private:
  

public:
  
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};
