#pragma once
#include "utils.hpp"
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
  public:
	simplexNode_P head = nullptr; //First simplex in the tree (0 vertex)
	simplexNode_P root = nullptr; //Empty node at root of tree (empty simplex)

	simplexTree(double, int);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexNode_P, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(simplexNode_P);
	void recurseInsert(simplexNode_P, unsigned, int, double, std::set<unsigned>);
	double findWeight(std::set<unsigned>);
	void deleteIndexRecurse(int, simplexNode_P);
	void deleteWeightEdgeGraph(int index);

	simplexNode_P find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, simplexNode_P);

	//virtual interface functions
	void outputComplex();
	double getSize();
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&);
	void deleteIterative(int);
	void deleteIndexRecurse(int);  // A wrapper for the actual deleteIndexRecurse method.
	void insert();
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	void prepareCofacets(int){return;}
	std::vector<simplexNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<simplexNode_P, simplexNode_P>&, bool);
	bool deletion(std::set<unsigned>);
	bool deletion(simplexNode_P);
	void expandDimensions(int){return;};
	void reduceComplex();
	void clear();
};

