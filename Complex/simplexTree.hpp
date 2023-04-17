#pragma once
#include "utils.hpp"
#include "simplexBase.hpp"
#include "kdTree.hpp"

#include <set>
#include <unordered_map>

// Header file for simplexTree class - see simplexTree.cpp for descriptions


template <class nodeType>
class simplexTree : public simplexBase<nodeType> {
  private:
  
  
  
  public:
  
    //SimplexTreeNode wraps a simplexNode and tree pointers; these are internal to the construction
    //  The child holds the shared pointer. Free children and siblings to clear the tree.
    template<typename nodetype>
    struct simplexTreeNode{

      struct cmpByIndex{
         bool operator()(const simplexTreeNode* lhs, const simplexTreeNode* rhs) const{
             return lhs->simpNode->index < rhs->simpNode->index;
         }
         bool operator()(const simplexTreeNode& lhs, const simplexTreeNode& rhs) const{
             return lhs.simpNode->index < rhs.simpNode->index;
         }
      };
      
      std::shared_ptr<nodeType> simpNode;                               //SharedPtr to the nodeType
      std::shared_ptr<simplexTreeNode<nodeType>> child = nullptr;       //SharedPtr to child
      std::shared_ptr<simplexTreeNode<nodeType>> sibling = nullptr;     //SharedPtr to sibling
      simplexTreeNode* parent = nullptr;
      std::set<simplexTreeNode*, cmpByIndex> children; 
      bool valid= true;	
          
      //SimplexNode Constructor Passthrough
      template<typename U = nodeType, typename std::enable_if<std::is_same<U, simplexNode>::value, int>::type = 0>
      simplexTreeNode(){simpNode = std::make_shared<nodeType>(simplexNode());}
      
      template<typename U = nodeType, typename std::enable_if<std::is_same<U, simplexNode>::value, int>::type = 0>
      simplexTreeNode(std::set<unsigned> simp, double wt){simpNode = std::make_shared<simplexNode>(simplexNode(simp, wt));}
      
      //AlphaNode Constructor Passthrough
      template<typename U = nodeType, typename std::enable_if<std::is_same<U, alphaNode>::value, int>::type = 0>
      simplexTreeNode(){simpNode = std::make_shared<alphaNode>(alphaNode());}
      
      template<typename U = nodeType, typename std::enable_if<std::is_same<U, alphaNode>::value, int>::type = 0>
      simplexTreeNode(std::set<unsigned> simp, double wt){simpNode = std::make_shared<alphaNode>(alphaNode(simp, wt));}
      
      //WitnessNode Constructor Passthrough
      template<typename U = nodeType, typename std::enable_if<std::is_same<U, witnessNode>::value, int>::type = 0>
      simplexTreeNode(){simpNode = std::make_shared<witnessNode>(witnessNode());}
      
      template<typename U = nodeType, typename std::enable_if<std::is_same<U, witnessNode>::value, int>::type = 0>
      simplexTreeNode(std::set<unsigned> simp, double wt){simpNode = std::make_shared<witnessNode>(witnessNode(simp, wt));}
      
    };
        
        
    
  
    //Typedefs for nodeType shared pointers; simplexTreeNode_P is the tree structure, templateNode_P is the simplex node structure
    typedef std::shared_ptr<simplexTreeNode<nodeType>> simplexTreeNode_P;
    typedef std::shared_ptr<nodeType> templateNode_P;
	simplexTreeNode_P root = nullptr; //Empty node at root of tree (empty simplex)
  
  
    //Presumably private functions (Not implemented in simplexBase)
	struct simplexTreeNode<nodeType>* find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, simplexTreeNode_P);
	struct simplexTreeNode<nodeType>* find(std::set<unsigned>::iterator, std::set<unsigned>::iterator, simplexTreeNode<nodeType>*);
  
  
	std::vector<std::set<templateNode_P, cmpByWeight<templateNode_P>>> getAllEdges();
	void recurseGetEdges(std::vector<std::set<templateNode_P, cmpByWeight<templateNode_P>>> &, simplexTreeNode_P, int, int);
  
    
	//Constructors
	simplexTree(double, int);
    
        
	std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexTreeNode_P, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);
	void printTree(simplexTreeNode_P);
	void printTree1(simplexTreeNode_P);
	void recurseInsert(simplexTreeNode<nodeType>*, unsigned, int, double, std::set<unsigned>);
	double findWeight(std::set<unsigned>);
	void deleteIndexRecurse(int, simplexTreeNode<nodeType>*);
	void deleteWeightEdgeGraph(int index);


	//virtual interface functions
    ~simplexTree(){};
	double getSize();
	int simplexCount();
	int vertexCount();
	void outputComplex();
    
    void insert();
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&);
    
	//bool find(std::vector<unsigned>);
	bool find(std::set<unsigned>);
    
	void deleteIterative(int);
	void deleteIndexRecurse(int);       // A wrapper for the actual deleteIndexRecurse method.
    
	void prepareCofacets(int){return;}
	void prepareFacets(int){return;}
    
	void expandDimensions(int){return;};
	std::vector<templateNode_P> expandDimension(std::vector<templateNode_P> edges){return std::vector<templateNode_P>();};
    
	//void reduceComplex();
    
    
    std::vector<templateNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<templateNode_P, templateNode_P>&, bool);
	std::vector<nodeType*> getAllCofacets(templateNode_P, const std::unordered_map<long long, templateNode_P>&, bool){return std::vector<nodeType*>();};
	std::vector<nodeType*> getAllCofacets(templateNode_P);
    std::vector<templateNode_P> getAllCofacets(const std::set<unsigned>&){return std::vector<templateNode_P>();};
    
	std::vector<templateNode_P> getAllDelaunayCofacets(templateNode_P){return std::vector<templateNode_P>();};
    
    std::vector<nodeType*> getAllFacets(nodeType*){return std::vector<nodeType*>();};
	std::vector<nodeType*> getAllFacets(templateNode_P);
	std::vector<templateNode_P> getAllFacets_P(templateNode_P);
	
    //Commenting out for now - not sure if the first argument type is what we want to use
    
	//std::vector<std::set<simplexNode_P, cmpByWeight>> getAllEdges();
	
    void recurseInsertDsimplex(simplexTreeNode_P node, std::vector<int> simp,std::vector<std::vector<double>> inputData);
    void buildAlphaComplex(std::vector<std::vector<int>> dsimplexmesh, int npts,std::vector<std::vector<double>> inputData);
	void validateNodes(simplexTreeNode_P headPointer);
	bool deletion(std::set<unsigned>);
	bool deletion(simplexTreeNode<nodeType>*);
	void graphInducedComplex(int dim,std::vector<std::vector<double>> inputData,double beta){return;};

	void reduceComplex();
	void clear();
    
    
    
};

