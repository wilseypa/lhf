#pragma once
#include <map>
#include <algorithm>
#include <set>
#include <iostream>
#include "utils.hpp"

// Header file for simplexBase class - see simplexTree.cpp for descriptions

class simplexBase {
  private:
  public:
	struct simplexNode{
		unsigned index;
		std::set<unsigned> simplex;
		simplexNode* child = nullptr;
		simplexNode* sibling = nullptr;
		simplexNode* parent = nullptr;
		double weight = 0;
	};
	std::vector<std::vector<simplexNode*>> simplexList;		//Holds ordered list of simplices in each dimension
																//Needs to sort by the weight for insertion
  
	long long nodeCount = 0;			//Total number of nodes stored
	long long indexCounter;				//Current insertion index
	
	utils ut;							//Utilities functions
	
	std::string simplexType = "simplexBase";
	double maxEpsilon;
	int maxDimension;
	std::vector<std::vector<double>> distMatrix;
	std::vector<std::vector<std::pair<std::vector<unsigned>, double>>> weightedGraph;
	int runningVectorCount = 0;
	std::vector<int> runningVectorIndices;
	int removedSimplices = 0;
	std::string stats = "RVIndex,Mean,Stdev,k,kNN_Mean,kNN_Stdev,Result\n";

	simplexBase();
	simplexBase(std::map<std::string, std::string>);
	simplexBase(double, int);
	void setConfig(std::map<std::string, std::string>);
	void setDistanceMatrix(std::vector<std::vector<double>> _distMatrix);
	simplexBase* newSimplex(const std::string &simplexT, std::map<std::string, std::string> configMap);
	bool (*streamEval) (std::vector<double>&, std::vector<std::vector<double>>&);
    bool streamEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
	void setStreamEvaluator(bool (*f) (std::vector<double>&, std::vector<std::vector<double>>&));

	//virtual interface functions
	virtual double getSize();
	virtual bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	virtual bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&, std::vector<double>&);
	virtual void deleteIterative(int);
	virtual void deleteIndexRecurse(int);  // A wrapper for the actual deleteIndexRecurse method.
	virtual void insert(std::vector<double>&);
	virtual bool find(std::vector<unsigned>);
	virtual bool find(std::set<unsigned>);
	virtual int simplexCount();
	virtual int vertexCount();
	virtual std::vector<simplexNode*> getDimEdges(int,double);
	virtual std::vector<std::vector<simplexNode*>> getAllEdges();
	virtual void expandDimensions(int);
	virtual void reduceComplex();
	virtual void clear();

	//Unused, possibly future
	virtual void outputSimplex();
};
