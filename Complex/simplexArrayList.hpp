#pragma once
#include "simplexBase.hpp"

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexArrayList : public simplexBase{
  private:
  public:
	simplexArrayList(double);
	std::string stats;	
	std::vector<std::vector<unsigned>> edges;
	std::vector<double> weights;
	
	double getSize();	
	std::vector<std::pair<double,std::vector<unsigned>>> getEdges(int,double);
	void insert(std::vector<double>);
	void find(std::vector<double>);
	int vertexCount();
	int simplexCount();
	void expandDimensions(int);
};

