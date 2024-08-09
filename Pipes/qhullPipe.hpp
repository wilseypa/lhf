#pragma once

#include <map>
#include "basePipe.hpp"

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"

template <typename nodeType>
class qhullPipe : public basePipe<nodeType>
{
private:
  std::string mode;

public:
  qhullPipe();
  void runPipe(pipePacket<nodeType> &inData);
  bool configPipe(std::map<std::string, std::string> &configMap);
  void outputData(pipePacket<nodeType> &);
  void qdelaunay_o(const orgQhull::Qhull &, std::vector<std::vector<unsigned>> &);
};
