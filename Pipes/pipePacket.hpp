#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for pipePacket class - see pipePacket.cpp for descriptions

struct bettiBoundaryTableEntry{
	int bettiDim;
	double birth;
	double death;
	std::set<unsigned> boundaryPoints;
}; 


class pipePacket {
  private:
  public:
	std::vector<bettiBoundaryTableEntry> bettiTable;
  
	pipePacket(const std::string &, const double, const int);
	pipePacket(const std::map<std::string, std::string>, const std::string&);
	std::string stats;
  
	std::vector<std::vector<double>> originalData;
	std::vector<unsigned> originalLabels;
	std::vector<std::vector<double>> fullData;
	simplexBase* complex;
	
	std::vector<std::vector<unsigned>> boundaries;
	std::set<double, std::greater<double>> weights;	
	std::string bettiOutput;
	
	double getSize();	
};

