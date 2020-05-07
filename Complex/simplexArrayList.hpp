#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

/*
class binomialTable{
	private:
		std::vector<std::vector<long long>> v; 
	public:
		binomialTable(unsigned n, unsigned k);
		long long binom(unsigned n, unsigned k);
};
*/

class simplexArrayList : public simplexBase{
  private:
	int indexCount;
	std::string stats;

	// long long ripsIndex(std::set<unsigned>& simplex, binomialTable& bin);
	// unsigned maxVertex(long long ripsIndex, unsigned high, unsigned low, unsigned k, binomialTable &bin);
	// std::vector<unsigned> getVertices(long long ripsIndex, int dim, unsigned n, binomialTable &bin);
  public:
	simplexArrayList(double, double, std::vector<std::vector<double>>*);
	double findWeight(std::set<unsigned>);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexNode*, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);

	//virtual interface functions
	double getSize();
	void insert(std::vector<double>&);
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	bool deletion(std::set<unsigned>);
	void expandDimensions(int);
	void reduceComplex();
	void clear();
};

