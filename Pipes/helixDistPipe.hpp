#pragma once
#include <map>
#include <set>
#include <vector>
#include "helixPipe.hpp"

// base constructor
template <typename nodeType>
class helixDistPipe : public helixPipe<nodeType>
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

public:
  helixDistPipe();
  void runPipe(pipePacket<nodeType> &inData);
  bool configPipe(std::map<std::string, std::string> &configMap);
  void outputData(pipePacket<nodeType> &);
};
