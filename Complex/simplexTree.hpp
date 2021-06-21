#pragma once
#include "utils.hpp"
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class simplexTree : public simplexBase {
  private:
	// For generation of combinations n choose r
	struct c_unique {
	unsigned current;
	c_unique() {current=0;}
	unsigned operator()() {return ++current;}
	} UniqueNumber;


	//SimplexTreeNode wraps a simplexNode and tree pointers
	struct simplexTreeNode{
		struct cmpByIndex{
			bool operator()(const simplexTreeNode* lhs, const simplexTreeNode* rhs) const{
				return lhs->simpNode->index < rhs->simpNode->index;
			}
			bool operator()(const simplexTreeNode& lhs, const simplexTreeNode& rhs) const{
				return lhs.simpNode->index < rhs.simpNode->index;
			}
		};
		
		simplexTreeNode* child = nullptr;
		simplexTreeNode* sibling = nullptr;
		simplexTreeNode* parent = nullptr;
		std::set<simplexTreeNode*, cmpByIndex> children;
		simplexNode_P simpNode;
		bool valid= true;	
		
		
		simplexTreeNode(){simpNode = std::make_shared<simplexNode>(simplexNode());}
		simplexTreeNode(std::set<unsigned> simp, double wt){simpNode = std::make_shared<simplexNode>(simplexNode(simp, wt));}
	};
  
	typedef std::shared_ptr<simplexTreeNode> simplexTreeNode_P;
  
	simplexTreeNode* find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, simplexTreeNode*);
	void recurseGetEdges(std::vector<std::set<simplexNode_P, cmpByWeight>> &, simplexTreeNode*, int, int);
  
  public:
	simplexTreeNode* root = nullptr; //Empty node at root of tree (empty simplex)

	simplexTree(double, int);
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexTreeNode*, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(simplexTreeNode*);
	void printTree1(simplexTreeNode*);
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
	void prepareFacets(int){return;}
	std::vector<simplexNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<simplexNode_P, simplexNode_P>&, bool = true);
	std::vector<simplexNode*> getAllCofacets(simplexNode_P, const std::unordered_map<long long, simplexNode_P>&, bool = true);
	std::vector<simplexNode*> getAllCofacets(simplexNode_P);
	std::vector<simplexNode_P> getAllDelaunayCofacets(simplexNode_P){return std::vector<simplexNode_P>();};
	std::vector<std::set<simplexNode_P, cmpByWeight>> getAllEdges();
	void recurseInsertDsimplex(simplexTreeNode* node, std::vector<unsigned> simp,std::vector<std::vector<double>> inputData);
    void buildAlphaComplex(std::vector<std::vector<int>> dsimplexmesh, int npts,std::vector<std::vector<double>> inputData);
	std::vector<simplexNode*> getAllFacets(simplexNode*);
	std::vector<simplexNode_P> getAllFacets_P(simplexNode_P);
	void validateNodes(simplexTreeNode* headPointer);
	bool deletion(std::set<unsigned>);
	bool deletion(simplexTreeNode*);
	void expandDimensions(int){return;};
	void checkInsertDsimplex(std::vector<unsigned> dsimplex,std::vector<std::vector<double>> inputData,double beta);
	void graphInducedComplex(std::vector<std::vector<double>> inputData,double beta);

	void reduceComplex();
	void clear();
};

