#pragma once
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
	bool isSorted = false;
	int nodeCount;
	std::vector<std::vector<std::pair<std::set<unsigned>, double>>> weightEdgeGraph;

	// void insertInductive();

  public:
	treeNode* head = nullptr;
	treeNode* root = nullptr;
	std::vector<treeNode*> dimensions;

	simplexTree(std::vector<std::vector<double>>);
	simplexTree(double, std::vector<std::vector<double>>, int);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(std::pair<std::set<unsigned>,double>, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(treeNode*);
	void recurseInsert(treeNode*, unsigned, int, double, std::set<unsigned>);
	double findWeight(std::set<unsigned>);
	void deleteIndexRecurse(int, treeNode*);
	void deleteWeightEdgeGraph(int index);
	treeNode* find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, treeNode*);

	//virtual interface functions
	double getSize();
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&, std::vector<double>&);
	void deleteIterative(int);
	void deleteIndexRecurse(int);  // A wrapper for the actual deleteIndexRecurse method.
	void insert(std::vector<double>&);
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	std::vector<std::vector<unsigned>> getDimEdges(int,double);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	std::vector<treeNode*> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<treeNode*, unsigned>&, bool);
	bool deletion(std::set<unsigned>);
	bool deletion(treeNode*);
	void expandDimensions(int){return;};
	void reduceComplex();
	std::vector<std::pair<double, std::vector<unsigned>>> getd0Pairs();
	void clear();
};

