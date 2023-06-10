#pragma once

#include "basePipe.hpp"
#include "utils.hpp"
#include "readInput.hpp"

// basePipe constructor
template <typename nodeType> 
class incrementalPipe : public basePipe<nodeType> {
  public:
    incrementalPipe();
    void runPipe(pipePacket<nodeType>& inData);
    bool configPipe(std::map<std::string, std::string> &configMap);
    void outputData(pipePacket<nodeType>&);
    // std::vector<std::vector<int>> qdelaunay_o(const delaunay &delaunay);
};
