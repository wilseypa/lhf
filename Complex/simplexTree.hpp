#pragma once
#include "utils.hpp"
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
	
	//SimplexTreeNode wraps a simplexNode and tree pointers
	struct simplexTreeNode{
		simplexTreeNode* child = nullptr;
		simplexTreeNode* sibling = nullptr;
		simplexTreeNode* parent = nullptr;
		simplexNode_P simpNode;
		
		
		simplexTreeNode(){simpNode = std::make_shared<simplexNode>(simplexNode());}
		simplexTreeNode(std::set<unsigned> simp, double wt){simpNode = std::make_shared<simplexNode>(simplexNode(simp, wt));}
	};
  
	typedef std::shared_ptr<simplexTreeNode> simplexTreeNode_P;
  
	simplexTreeNode* find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, simplexTreeNode*);
  
  public:
	simplexTreeNode* root = nullptr; //Empty node at root of tree (empty simplex)

	simplexTree(double, int);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexTreeNode*, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(simplexTreeNode*);
	void recurseInsert(simplexTreeNode*, unsigned, int, double, std::set<unsigned>);
	double findWeight(std::set<unsigned>);
	void deleteIndexRecurse(int, simplexTreeNode*);
	void deleteWeightEdgeGraph(int index);


	//virtual interface functions
	void outputComplex();
	double getSize();
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&);
	void deleteIterative(int);
	void deleteIndexRecurse(int);  // A wrapper for the actual deleteIndexRecurse method.
	void insert();
	bool find(std::set<unsigned>);
	int simplexCount();
	int vertexCount();
	void prepareCofacets(int){return;}
	std::vector<simplexNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<simplexNode_P, simplexNode_P>&, bool);
	bool deletion(std::set<unsigned>);
	bool deletion(simplexTreeNode*);
	void expandDimensions(int){return;};
	void reduceComplex();
	void clear();
};

