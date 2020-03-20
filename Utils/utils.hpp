#pragma once

#ifndef UTILS_HPP_INCL
#define UTILS_HPP_INCL

#include <set>

// Header file for utils class - see utils.cpp for descriptions

class utils {
  private:
	std::string debug = "0";
	std::string outputFile = "console";
  
  public:
	utils();
	utils(std::string, std::string);
	double computeMaxRadius(int, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<unsigned>);
	std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> separatePartitions(int, std::vector<std::vector<double>>, std::vector<unsigned>);
	std::vector<std::vector<std::vector<double>>> separateBoundaryPartitions(std::vector<std::set<unsigned>>, std::vector<std::vector<double>>, std::vector<unsigned>);
	std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> separatePartitions(double, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<unsigned>);
	void print2DVector(const std::vector<std::vector<unsigned>>&);
	void print1DVector(const std::vector<unsigned>&);
	void print1DVector(const std::set<unsigned>&);
	void print1DVector(const std::vector<double>&);
	std::vector<double> feature_distance(std::vector<double>*, std::vector<double>*);
	double vectors_distance(const double&, const double&);
	double vectors_distance(const std::vector<double>&, const std::vector<double>&);
	void print1DSet(const auto&);	
	std::set<unsigned> setXOR(std::set<unsigned>&, std::set<unsigned>&);
	std::set<unsigned> setIntersect(std::set<unsigned>, std::set<unsigned>, bool isSorted);
	std::vector<unsigned> setIntersect(std::vector<unsigned>, std::vector<unsigned>, bool);
	std::vector<std::set<unsigned>> getSubsets(std::set<unsigned>, int);
	std::vector<unsigned> symmetricDiff(std::vector<unsigned>, std::vector<unsigned>, bool);
	std::vector<unsigned> symmetricDiff(std::set<unsigned>, std::set<unsigned>, bool);
	std::vector<unsigned> setUnion(std::vector<unsigned>, std::vector<unsigned>, bool);
	std::pair<std::vector<unsigned>, std::vector<unsigned>> intersect(std::vector<unsigned>, std::vector<unsigned>, bool);
	
	//Utility functions for writing to console/debug file
	void writeLog(std::string module, std::string message);
	void writeDebug(std::string module, std::string message);
	void writeError(std::string module, std::string error){writeLog(module,error);return;};
	void writeFile(std::string fullMessage);
	
	static bool sortBySecond(const std::pair<std::set<unsigned>, double> &, const std::pair<std::set<unsigned>, double> &);
	std::vector<std::set<unsigned>> getSubsets(std::set<unsigned> set);
	std::vector<std::vector<unsigned>> getSubsets(std::vector<unsigned> set);
	
	std::vector<double> nearestNeighbors(std::vector<double>&, std::vector<std::vector<double>>&);
};

#endif
