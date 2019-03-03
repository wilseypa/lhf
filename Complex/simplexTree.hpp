#pragma once
#include "simplexBase.hpp"

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
	int maxDim;
	struct treeNode{
		int index;
		treeNode* child = nullptr;
		treeNode* sibling = nullptr;
		treeNode* parent = nullptr;
	};
			
	std::vector<treeNode*> dimensions;		
	
	int indexCounter;
	treeNode* head;
	int nodeCount;
  
  public:
	simplexTree(std::vector<std::vector<double>>);
	simplexTree(double, std::vector<std::vector<double>>, int);
	
	bool isLeaf;
	
	double getSize();
	void printTree(treeNode*);

	// At this time, let's just assume that each simplex is labeled by a key that
	// can, in general, be considered as a string.
	void insert(std::vector<double>);
	bool deletion(simplexTree*&, std::string);
	bool search(std::string);
	bool haveChild(simplexTree const*);
	int vertexCount();
	int simplexCount();
};

