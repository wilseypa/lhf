#pragma once
#include "simplexBase.hpp"

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
  public:
	simplexTree();
	
	double getSize();	
	void insert(std::vector<double>);
	void find(std::vector<double>);
	int vertexCount();
	int simplexCount();
};

