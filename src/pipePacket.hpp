#pragma once
#include "simplexBase.hpp"

// Header file for pipePacket class - see pipePacket.cpp for descriptions

class pipePacket {
  private:
  public:
	pipePacket(const std::string &);
	std::string stats;
  
	struct pipeData{
		std::vector<std::vector<double>> originalData;
		std::vector<std::vector<double>> workingData;
		simplexBase* complex;
		std::vector<std::vector<unsigned>> edges;
		std::vector<double> weights;
	} workData;
	
	
	double getSize();	
};

