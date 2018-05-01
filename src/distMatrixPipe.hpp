#pragma once

// Header file for distMatrixPipe class - see distMatrixPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"

class distMatrixPipe : public basePipe {
  private:
  public:
    distMatrixPipe();
    std::vector<std::vector<double>> runPipe(std::vector<std::vector<double>> inData);
    bool configPipe(std::map<std::string, std::string> configMap);
};

