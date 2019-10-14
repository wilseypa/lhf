#pragma once

// Header file for kdTree class - see kdTree.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

class kdTree : public preprocessor
{
private:
  

public:
  
    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};
