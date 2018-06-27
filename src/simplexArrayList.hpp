#pragma once
#include "simplexBase.hpp"

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexArrayList : public simplexBase{
  private:
  public:
	simplexArrayList();
	std::string stats;	
	std::vector<std::vector<unsigned>> edges;
	std::vector<double> weights;
	
	double getSize();	
	void insert(std::vector<double>);
	void find(std::vector<double>);
	int vertexCount();
	int simplexCount();
};

