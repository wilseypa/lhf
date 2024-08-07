#pragma once

#include "basePipe.hpp"
#include "utils.hpp"

// basePipe constructor
template <typename nodeType>
class incrementalPipe : public basePipe<nodeType>
{
private:
  std::vector<std::vector<double>> &inputData;
  std::map<std::vector<short>, short> inner_d_1_shell;
  std::set<std::vector<short>> dsimplexes;
  std::vector<short> search_space;
  unsigned dim;
  unsigned data_set_size;
  std::vector<double> solvePlaneEquation(const std::vector<short> &points);
  std::vector<short> first_simplex();
  void cospherical_handler(std::vector<short> &simp, int &tp, short &omission, std::vector<std::vector<double>> &distMatrix);
  int expand_d_minus_1_simplex(std::vector<short> &simp_vector, short &omission, std::vector<std::vector<double>> &distMatrix);

public:
  incrementalPipe();
  void runPipe(pipePacket<nodeType> &inData);
  bool configPipe(std::map<std::string, std::string> &configMap);
  void outputData(pipePacket<nodeType> &);
};
