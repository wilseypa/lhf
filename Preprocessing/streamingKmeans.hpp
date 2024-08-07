#pragma once

// Header file for streamingkMeans class - see streamingkMeans.cpp for descriptions
#include <map>
#include "preprocessor.hpp"

template <typename nodeType>
class streamingKmeans : public preprocessor<nodeType>
{
private:
	int numClusters;
	int num_iterations;

public:
	streamingKmeans();
	double dotProd(const std::vector<double> &a, const std::vector<double> &b);
	double dotProd2D(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &b);
	std::vector<double> approxNearestNeighbor(std::vector<std::vector<double>> &facilities, std::vector<std::pair<double, int>> &sortedApproxFacils, std::vector<double> omega, int x, int size, pipePacket<nodeType>(inData));
	std::vector<double> approxHat(std::vector<std::vector<double>> &summedCentroidVectors, std::vector<std::pair<double, int>> &sortedApproxFacilsHat, std::vector<double> omega, int xHat, int size);
	int binarySearch(std::vector<std::pair<double, int>> &sorted, std::vector<double> omega, int n, double target);
	double randDouble(); // Returns a random double in the range [0,1)
	bool prob(double f);
	int random(int low, int high); // returns random int in the range [low, high)

	void runPreprocessor(pipePacket<nodeType> &);
	bool configPreprocessor(std::map<std::string, std::string> &);
};
