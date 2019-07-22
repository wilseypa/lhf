#pragma once
#include "simplexBase.hpp"
#include <set>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
	bool isSorted = false;
	
	struct cmp{
		bool operator()(const std::pair<std::set<unsigned>, double> &a, const std::pair<std::set<unsigned>, double> &b){
			return (a.second < b.second);
		}
	};
  
	int maxDim;
	struct treeNode{
		unsigned index;
		treeNode* child = nullptr;
		treeNode* sibling = nullptr;
		treeNode* parent = nullptr;
		double weight = 0;
	};
	
	std::vector<std::vector<std::pair<std::set<unsigned>, double>>> weightEdgeGraph;
			
	std::vector<treeNode*> dimensions;		
	
	unsigned indexCounter;
	treeNode* head;
	int nodeCount;
  
  public:
	simplexTree(std::vector<std::vector<double>>);
	simplexTree(double, std::vector<std::vector<double>>, int);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	
	bool isLeaf;
	
	double getSize();
	void printTree(treeNode*);

	// At this time, let's just assume that each simplex is labeled by a key that
	// can, in general, be considered as a string.
	void insert(std::vector<double>&);
	bool deletion(simplexTree*&, std::string);
	bool search(std::string);
	bool haveChild(simplexTree const*);
	void insertInductive();
	void recurseInsert(treeNode*, unsigned, int, double, std::set<unsigned>);
	int vertexCount();
	int simplexCount();
	void reduceComplex();
};

