#pragma once
#include <set>
#include "simplexBase.hpp"

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class indSimplexTree : public simplexBase {
  private:
	bool isSorted = false;
	
	struct cmp{
		bool operator()(const std::pair<std::set<unsigned>, double> &a, const std::pair<std::set<unsigned>, double> &b){
			return (a.second < b.second);
		}
	};
  
	int maxDim;
	
	std::vector<std::vector<indTreeNode*>> dimensions;		
	
	unsigned indexCounter;
	indTreeNode* head;
	int nodeCount;
  
  public:
	std::vector<std::vector<graphEntry>> indexedGraph;
	bool isLeaf;
	
	indSimplexTree(std::vector<std::vector<double>>);
	indSimplexTree(double, std::vector<std::vector<double>>, int);
	std::set<unsigned> getFaces(graphEntry ge);
	std::vector<std::vector<graphEntry>> coreduction(graphEntry);
	static bool compareByWeight(const graphEntry &, const graphEntry &);
	void printTree(indTreeNode*);
	bool search(std::set<unsigned>);
	bool haveChild(indSimplexTree const*);
	double getWeight(std::set<unsigned>);
	void sortAndBuildGraph();
	void insertInductive();
	void recurseInsert(indTreeNode*, unsigned, int, double, std::set<unsigned>);
	std::pair<std::vector<std::set<unsigned>>,std::vector<std::set<unsigned>>> recurseReduce(std::set<unsigned>,int,std::set<unsigned>,int, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	
	//virtual interface functions
	double getSize();
	void insert(std::vector<double>&);
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	std::vector<std::vector<unsigned>> getDimEdges(int,double);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	std::vector<std::vector<graphEntry>> getIndexEdges(double);
	bool deletion(indTreeNode*);
	bool deletion(std::set<unsigned>);
	void expandDimensions(int){return;};
	void reduceComplex();
	
	
};

