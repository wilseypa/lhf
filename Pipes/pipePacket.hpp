#pragma once
#include "simplexBase.hpp"

// Header file for pipePacket class - see pipePacket.cpp for descriptions

class pipePacket {
  private:
  public:
	pipePacket(const std::string &, const double);
	std::string stats;
  
	struct pipeData{
		std::vector<std::vector<double>> originalData;
		std::vector<unsigned> originalLabels;
		std::vector<std::vector<double>> upscaleData;
		simplexBase* complex;
	} workData;
	
	std::vector<std::vector<unsigned>> boundaries;
	std::string bettiOutput;
	
	double getSize();	
};

