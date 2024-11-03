#pragma once

#include "basePipe.hpp"
#include "utils.hpp"

// basePipe constructor
template <typename nodeType>
class helixPipe : public basePipe<nodeType>
{
protected:
  std::vector<std::vector<double>> inputData;
  std::vector<std::vector<short>> dsimplexmesh;
  std::set<std::vector<short>> spherical_dsimplexes;
  std::vector<short> search_space;
  unsigned dim;
  unsigned data_set_size;
  std::vector<double> solvePlaneEquation(const std::vector<short> &points);
  std::vector<short> first_simplex();
  void reduce(std::set<std::vector<short>> &outer_dsimplexes, std::vector<std::pair<std::vector<short>, short>> &inner_d_1_shell, std::vector<std::vector<short>> &dsimplexes);
  void cospherical_handler(std::vector<short> &simp, int &tp, short &omission, std::vector<std::vector<double>> &distMatrix);
  short expand_d_minus_1_simplex(std::vector<short> &simp_vector, short &omission, std::vector<std::vector<double>> &distMatrix);

public:
  helixPipe();
  void runPipe(pipePacket<nodeType> &inData);
  bool configPipe(std::map<std::string, std::string> &configMap);
  void outputData(pipePacket<nodeType> &);
};
