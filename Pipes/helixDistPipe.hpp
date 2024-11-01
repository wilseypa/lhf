#pragma once
#include <map>
#include <set>
#include <vector>
#include "basePipe.hpp"

// base constructor
template <typename nodeType>
class helixDistPipe : public basePipe<nodeType>
{
private:
  std::vector<std::vector<double>> inputData;
  std::vector<std::vector<double>> distMatrix;
  std::map<std::vector<short>, short> inner_d_1_shell_map;
  std::set<std::vector<short>> dsimplexes;
  std::set<std::vector<short>> spherical_dsimplexes;
  std::set<short> search_space;
  unsigned dim;
  unsigned data_set_size;
  std::vector<double> solvePlaneEquation(const std::vector<short> &points);
  std::vector<short> first_simplex();
  void cospherical_handler(std::vector<short> &simp, int &tp, short &omission);
  short expand_d_minus_1_simplex(std::vector<short> &simp_vector, short &omission);

public:
  helixDistPipe();
  void runPipe(pipePacket<nodeType> &inData);
};
