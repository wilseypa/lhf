#pragma once
#include "simplexBase.hpp"
#include <set>
#include "utils.hpp"

// Header file for pipePacket class - see pipePacket.cpp for descriptions



template<typename T>
class pipePacket {
  private:
  public:
	std::vector<bettiBoundaryTableEntry> bettiTable;
	std::string ident;
  
	pipePacket<T>(const std::string &, const double, const int);
	pipePacket<T>(std::map<std::string, std::string>, const std::string&);
	std::string stats;
	std::string runLog;
  
	std::vector<std::vector<double>> workData;
	std::vector<unsigned> centroidLabels;
	std::vector<std::vector<double>> inputData;
	std::vector<std::vector<double>> distMatrix;
	simplexBase<T>* complex = nullptr;
	
	std::vector<std::set<unsigned>> boundaries;
	std::set<double, std::greater<double>> weights;	
	std::string bettiOutput;
	
	double getSize();	
	std::string getStats();
};


//Explicit Template Class Instantiation
template class pipePacket<simplexNode>;
template class pipePacket<alphaNode>;
