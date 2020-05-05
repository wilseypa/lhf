#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
	bool isSorted = false;

  public:

	simplexBase::simplexNode* head = nullptr;
	//simplexBase::simplexTree(std::vector<std::vector<double>>);

	simplexTree(double, std::vector<std::vector<double>>*, int);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexBase::simplexNode*, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(simplexBase::simplexNode*);
	void recurseInsert(simplexBase::simplexNode*, unsigned, int, double, std::set<unsigned>);
	double findWeight(std::set<unsigned>);
	void deleteIndexRecurse(int, simplexBase::simplexNode*); 
	void deleteWeightEdgeGraph(int index);
	simplexNode* find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, simplexNode*);
	std::vector<simplexNode*> getAllCofacets(const std::set<unsigned>&);

	//virtual interface functions
	double getSize();
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&, std::vector<double>&);
	void deleteIterative(simplexBase::simplexNode*);
	void deleteIndexRecurse(int);  // A wrapper for the actual deleteIndexRecurse method.
	void insert(std::vector<double>&);
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	bool deletion(std::set<unsigned>);
	bool deletion(simplexBase::simplexNode*);
	void expandDimensions(int){return;};
	void reduceComplex();
	void clear();
};

