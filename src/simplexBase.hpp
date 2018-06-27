#pragma once

// Header file for simplexBase class - see simplexTree.cpp for descriptions

class simplexBase {
  private:
  public:
	std::string simplexType;
	simplexBase();
	simplexBase* newSimplex(const std::string &simplexT);
	
	double getSize();
	void insert(std::vector<double>);
	void find(std::vector<double>);
	int simplexCount();
	int vertexCount();
};

