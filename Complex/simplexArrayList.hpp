#pragma once
#include "simplexBase.hpp"

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexArrayList : public simplexBase{
  private:
	int indexCount;
  public:
	simplexArrayList(double, std::vector<std::vector<double>>);
	std::string stats;
	
	double getSize();	
	std::vector<std::vector<unsigned>> getEdges(int,double);
	std::vector<std::vector<std::vector<unsigned>>> getAllEdges(double);
	void insert(std::vector<double>);
	void find(std::vector<double>);
	int vertexCount();
	int simplexCount();
	void expandDimensions(int);
};

