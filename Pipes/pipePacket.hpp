#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for pipePacket class - see pipePacket.cpp for descriptions

class pipePacket {
  private:
  public:
	pipePacket(const std::string &, const double, const int);
	pipePacket(const std::map<std::string, std::string>, const std::string&);
	std::string stats;
  
	std::vector<std::vector<double>> originalData;
	std::vector<unsigned> originalLabels;
	std::vector<std::vector<double>> reducedData;
	simplexBase* complex;
	
	std::vector<std::vector<unsigned>> boundaries;
	std::set<double, std::greater<double>> weights;	
	std::string bettiOutput;
	
	double getSize();	
};

