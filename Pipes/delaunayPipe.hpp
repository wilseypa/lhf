#pragma once

#include "basePipe.hpp"

#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

// basePipe constructor
template <typename nodeType>
class delaunayPipe : public basePipe<nodeType>
{
private:
public:
  delaunayPipe();
  void runPipe(pipePacket<nodeType> &inData);
  bool configPipe(std::map<std::string, std::string> &configMap);
  void outputData(pipePacket<nodeType> &);
  // std::vector<std::vector<int>> qdelaunay_o(const delaunay &delaunay);
};
