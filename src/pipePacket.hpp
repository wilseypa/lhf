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
		simplexBase* complex;
	} workData;
	
	
	double getSize();	
};

