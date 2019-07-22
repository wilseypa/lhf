#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexArrayList : public simplexBase{
  private:
	int indexCount;
	std::string stats;
  public:
	simplexArrayList(double, std::vector<std::vector<double>>);
		
	//virtual interface functions
	double getSize();
	void insert(std::vector<double>&);
	void find(std::vector<double>);
	int simplexCount();
	int vertexCount();
	std::vector<std::vector<unsigned>> getDimEdges(int,double);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	bool deletion(std::vector<unsigned>);
	void expandDimensions(int);
	void reduceComplex();
};

