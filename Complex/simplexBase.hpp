#pragma once
#include <set>

// Header file for simplexBase class - see simplexTree.cpp for descriptions

class simplexBase {
  private:
  public:
	std::string simplexType;
	simplexBase();
	simplexBase(double, int);
	void setDistanceMatrix(std::vector<std::vector<double>> _distMatrix);
	simplexBase* newSimplex(const std::string &simplexT);
	double maxEpsilon;
	int maxDimension;
	std::vector<std::vector<double>> distMatrix;
	std::vector<std::vector<std::vector<unsigned>>> weightedGraph;
	
	virtual double getSize();
	virtual void insert(std::vector<double>&);
	virtual void find(std::vector<double>);
	virtual int simplexCount();
	virtual int vertexCount();
	virtual std::vector<std::vector<unsigned>> getEdges(int,double);
	virtual std::vector<std::vector<std::pair<std::set<unsigned>, double>>> getAllEdges(double);
	virtual void outputSimplex();
	virtual void expandDimensions(int);
};
