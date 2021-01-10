#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for pipePacket class - see pipePacket.cpp for descriptions




class pipePacket {
  private:
  public:
  	std::set<std::vector<std::set<unsigned>>> VUdsimplexes;  // Keep all the unique and valid d-simplexes enumerations
    std::set<bettiTables,cmpByPandE> bTbs;	//Store betties from individual VSC

	std::vector<bettiBoundaryTableEntry> bettiTable;
  
	pipePacket(const std::string &, const double, const int);
	pipePacket(std::map<std::string, std::string>, const std::string&);
	std::string stats;
  
	std::vector<std::vector<double>> workData;
	std::vector<unsigned> centroidLabels;
	std::vector<std::vector<double>> inputData;
	std::vector<std::vector<double>> distMatrix;
	simplexBase* complex = nullptr;
	
	std::vector<std::set<unsigned>> boundaries;
	std::set<double, std::greater<double>> weights;	
	std::string bettiOutput;
	
	double getSize();	
};

