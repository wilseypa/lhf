#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
	bool isSorted = false;
	unsigned indexCounter;
	int nodeCount;
	std::vector<std::vector<std::pair<std::set<unsigned>, double>>> weightEdgeGraph;

	void insertInductive();

  public:
	struct treeNode{
		unsigned index;
		std::set<unsigned> simplex;
		treeNode* child = nullptr;
		treeNode* sibling = nullptr;
		treeNode* parent = nullptr;
		double weight = 0;
	};

	treeNode* head = nullptr;
	std::vector<treeNode*> dimensions;

	simplexTree(std::vector<std::vector<double>>);
	simplexTree(double, std::vector<std::vector<double>>, int);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(std::pair<std::set<unsigned>,double>, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(treeNode*);
	void recurseInsert(treeNode*, unsigned, int, double, std::set<unsigned>);
	double findWeight(std::set<unsigned>);
	void deleteIndexRecurse(int, treeNode*);

	//virtual interface functions
	double getSize();
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	void deleteIterative(int, int);
	void insert(std::vector<double>&);
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	std::vector<std::vector<unsigned>> getDimEdges(int,double);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	bool deletion(std::set<unsigned>);
	bool deletion(treeNode*);
	void expandDimensions(int){return;};
	void reduceComplex();
	std::vector<std::pair<double, std::vector<unsigned>>> getd0Pairs();
	void clear();
};

