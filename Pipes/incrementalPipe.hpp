#pragma once

#include "basePipe.hpp"
#include "utils.hpp"
#include "readInput.hpp"

// basePipe constructor
template <typename nodeType> 
class incrementalPipe : public basePipe<nodeType> {
  private:
  	std::vector<std::vector<double>>& inputData;
    std::map<std::vector<short>, short> inner_d_1_shell;
    std::vector<short> search_space;
    unsigned dim;
    unsigned data_set_size;
  public:
    incrementalPipe();
    void runPipe(pipePacket<nodeType>& inData);
    std::vector<double> solvePlaneEquation(const std::vector<short> &points);
    short validate(std::vector<short>& simp, short &triangulation_point, short &omission);
    std::vector<short> first_simplex();
    int expand_d_minus_1_simplex(std::vector<short> &simp_vector, short &omission, std::vector<std::vector<double>>& distMatrix);
    bool configPipe(std::map<std::string, std::string> &configMap);
    void outputData(pipePacket<nodeType>&);
    // std::vector<std::vector<int>> qdelaunay_o(const delaunay &delaunay);
};
