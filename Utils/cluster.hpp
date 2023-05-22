#pragma once

#include <set>
#include <vector>
#include "utils.hpp"

class kmeansplusplus {
  private:
    utils ut;

  public:
	kmeansplusplus(){ ut = utils();};

	void clusterData(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<unsigned> &, int, int, int);



};
