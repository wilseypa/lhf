#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexArrayList : public simplexBase{
  private:
	int indexCount;
	std::string stats;
  public:
	simplexArrayList(double, double, std::vector<std::vector<double>>);
	double findWeight(std::vector<unsigned>);
	std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<unsigned>>> recurseReduce(std::pair<std::vector<unsigned>,double>, std::vector<std::vector<unsigned>>, std::vector<std::vector<unsigned>>);

		
	//virtual interface functions
	double getSize();
	void insert(std::vector<double>&);
	bool find(std::vector<unsigned>);
	int simplexCount();
	int vertexCount();
	std::vector<std::vector<unsigned>> getDimEdges(int,double);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	bool deletion(std::vector<unsigned>);
	void expandDimensions(int);
	void reduceComplex();
	void clear();
};

