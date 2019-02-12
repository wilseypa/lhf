#pragma once

// Header file for simplexBase class - see simplexTree.cpp for descriptions

class simplexBase {
  private:
  public:
	std::string simplexType;
	simplexBase();
	simplexBase(double);
	void setDistanceMatrix(std::vector<std::vector<double>> _distMatrix);
	simplexBase* newSimplex(const std::string &simplexT);
	double maxEpsilon;
	std::vector<std::vector<double>> distMatrix;
	std::vector<std::vector<std::pair<double,std::vector<unsigned>>>> weightedGraph;
	
	std::vector<double> weights;	
	
	virtual double getSize();
	virtual void insert(std::vector<double>);
	virtual void find(std::vector<double>);
	virtual int simplexCount();
	virtual int vertexCount();
	virtual std::vector<std::vector<unsigned>> getEdges(int,double);
	virtual std::vector<std::pair<double,std::vector<unsigned>>> getAllEdges();
	virtual void outputSimplex();
	virtual void expandDimensions(int);
};
