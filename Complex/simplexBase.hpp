#pragma once

// Header file for simplexBase class - see simplexTree.cpp for descriptions

class simplexBase {
  private:
  public:
	std::string simplexType;
	simplexBase();
	simplexBase(double);
	simplexBase* newSimplex(const std::string &simplexT);
	double maxEpsilon;
	std::vector<std::vector<double>> data;
	std::vector<std::pair<double,std::vector<unsigned>>> weightedGraph;
	
	//Temporary!:
	std::vector<double> weights;
	std::vector<std::vector<unsigned>> edges;
	
	
	virtual double getSize();
	virtual void insert(std::vector<double>);
	virtual void find(std::vector<double>);
	virtual int simplexCount();
	virtual int vertexCount();
	virtual std::vector<std::vector<unsigned>> getEdges(int,double);
	virtual void outputSimplex();
	virtual void expandDimensions(int);
};

