#pragma once
#include "simplexBase.hpp"

// Define the maximum number of points whose corresponding simplices can be stored
// in the simplex tree.
#define MAX_POINTS 2000

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
  public:
	simplexTree();
	bool isLeaf;
	simplexTree* character[MAX_POINTS];
	
	double getSize();

	// At this time, let's just assume that each simplex is labeled by a key that
	// can, in general, be considered as a string.
	void insert(std::string);
	bool deletion(simplexTree*&, std::string);
	bool search(std::string);
	bool haveChild(simplexTree const*);
	int vertexCount();
	int simplexCount();
};

