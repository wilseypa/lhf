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
  
	indSimplexTree(std::vector<std::vector<double>>);
	indSimplexTree(double, std::vector<std::vector<double>>, int);
	std::vector<std::vector<graphEntry>> getIndexEdges(double);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	std::set<unsigned> getFaces(graphEntry ge);
	
	static bool compareByWeight(const graphEntry &, const graphEntry &);
	bool isLeaf;
	
	double getSize();
	void printTree(indTreeNode*);

	// At this time, let's just assume that each simplex is labeled by a key that
	// can, in general, be considered as a string.
	void insert(std::vector<double>&);
	bool deletion(indSimplexTree*&, std::string);
	bool search(std::set<unsigned>);
	unsigned find(std::set<unsigned>);
	bool haveChild(indSimplexTree const*);
	void insertInductive();
	void recurseInsert(indTreeNode*, unsigned, int, double, std::set<unsigned>);
	int vertexCount();
	int simplexCount();
	void sortAndBuildGraph();
};

