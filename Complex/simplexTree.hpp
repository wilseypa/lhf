#pragma once
#include "utils.hpp"
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
  public:
	simplexNode* head = nullptr; //First simplex in the tree (0 vertex)
	simplexNode* root = nullptr; //Empty node at root of tree (empty simplex)

	simplexTree(double, std::vector<std::vector<double>>*, int);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexNode*, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(simplexNode*);
	void recurseInsert(simplexNode*, unsigned, int, double, std::set<unsigned>);
	double findWeight(std::set<unsigned>);
	void deleteIndexRecurse(int, simplexNode*); 
	void deleteWeightEdgeGraph(int index);

	simplexNode* find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, simplexNode*);

	//virtual interface functions
	void outputComplex();
	double getSize();
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&, std::vector<double>&);
	void deleteIterative(simplexNode*);
	void deleteIndexRecurse(int);  // A wrapper for the actual deleteIndexRecurse method.
	void insert();
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	void prepareCofacets(int){return;}
	std::vector<simplexNode*> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<simplexNode*, simplexNode*>&, bool);
	bool deletion(std::set<unsigned>);
	bool deletion(simplexNode*);
	void expandDimensions(int){return;};
	void reduceComplex();
	void clear();
};

