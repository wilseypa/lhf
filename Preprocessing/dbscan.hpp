#pragma once

// Header file for dbscan class - see dbscan.cpp for descriptions
#include <map>
#include <vector>
#include "kdTree.hpp"

class dbscan
{
private:
  static void expandCluster(const std::vector<std::vector<double>> &data,
                            std::vector<int> &labels,
                            std::vector<size_t> &neighbors,
                            int clusterLabel,
                            kdTree &tree,
                            int minPoints,
                            double epsilon);

public:
  static std::vector<int> cluster(const std::vector<std::vector<double>> &data, int minPoints, double epsilon, int size);
};
