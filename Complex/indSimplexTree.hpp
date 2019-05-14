#pragma once
#include "simplexBase.hpp"
#include <set>

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
	struct indTreeNode{
		unsigned index;
		indTreeNode* child = nullptr;
		indTreeNode* sibling = nullptr;
		indTreeNode* parent = nullptr;
		std::set<unsigned> simplexSet;
		double weight = 0;
	};
	
	struct graphEntry{
		std::set<unsigned> simplexSet;
		double weight = 0;
		indTreeNode* entry = nullptr;	
		
		graphEntry(){}
		graphEntry(std::set<unsigned> simp, double wt, indTreeNode* ent) { 
			simplexSet = simp; weight = wt; entry = ent;
		}
	};
	
	std::vector<std::vector<graphEntry>> indexedGraph;
	std::vector<int> dimCounts = {6, 10, 20};
			
	std::vector<indTreeNode*> dimensions;		
	
	unsigned indexCounter;
	indTreeNode* head;
	int nodeCount;
  
  public:
	indSimplexTree(std::vector<std::vector<double>>);
	indSimplexTree(double, std::vector<std::vector<double>>, int);
	std::vector<std::vector<indSimplexTree::graphEntry>> getIndexEdges(double);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> getAllEdges(double);
	
	bool isLeaf;
	
	double getSize();
	void printTree(indTreeNode*);

	// At this time, let's just assume that each simplex is labeled by a key that
	// can, in general, be considered as a string.
	void insert(std::vector<double>&);
	bool deletion(indSimplexTree*&, std::string);
	bool search(std::set<unsigned>);
	bool haveChild(indSimplexTree const*);
	void insertInductive();
	void recurseInsert(indTreeNode*, unsigned, int, double, std::set<unsigned>);
	int vertexCount();
	int simplexCount();
	void sortAndBuildGraph();
};

