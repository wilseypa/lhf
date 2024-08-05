#pragma once
#include <map>
#include <algorithm>
#include <set>
#include <memory>
#include <iostream>
#include <unordered_map>
#include "utils.hpp"

// Header file for simplexBase class - see simplexTree.cpp for descriptions

template <typename nodeType>
struct cmpByWeight{
	bool operator()(nodeType a, nodeType b) const{
		if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
			auto itA = a->simplex.rbegin(), itB = b->simplex.rbegin();
			while(itA != a->simplex.rend()){
				if(*itA != *itB) return *itA > *itB;
				++itA; ++itB;
			}
			return false;
		} else{
			return a->weight < b->weight;
		}
	}
};


template <class nodeType>
class simplexBase {
  private:
  
  public:
	typedef std::shared_ptr<nodeType> templateNode_P;
	std::vector<std::set<templateNode_P, cmpByWeight<templateNode_P>>> simplexList;		//Holds ordered list of simplices in each dimension
	std::vector<std::vector<unsigned>> dsimplexmesh;
	unsigned simplexOffset = 0;

	
	long long nodeCount = 0;					//Total number of nodes stored
	long long indexCounter;						//Current insertion index

	utils ut;									//Utilities functions
	std::string simplexType = "simplexBase";	//Complex Type Identifier
	std::string complexType = "";
	std::string simplicialComplex = "";
	templateNode_P root;							//Root of the simplexNode tree (if applicable)
	templateNode_P head;							//Root of the simplexNode tree (if applicable)

	double maxEpsilon;							//Maximum epsilon, loaded from configuration
	int maxDimension;							//Maximum dimension, loaded from configuration
	double alphaFilterationValue;						//alpha FilterationValue for alpha Complex
	std::vector<std::vector<double>>* distMatrix;	//Pointer to distance matrix for current complex
    std::vector<std::vector<bool>>* incidenceMatrix;

	//For sliding window implementation, tracks the current vectors inserted into the window
	//		Note - these point to the d0 simplexNodes; index, weight, etc. can be obtained
	std::vector<int> runningVectorIndices;   // std::vector<simplexNode*> runningVectorIndices;
	int runningVectorCount = 0;		//How many total points have been inserted into complex?
									//		Is this different than indexCounter?

	int removedSimplices = 0;
	std::string stats = "RVIndex,Mean,Stdev,k,kNN_Mean,kNN_Stdev,Result\n";

	//Constructors
	simplexBase();
	simplexBase(std::map<std::string, std::string>&);
	simplexBase(double, int);

	//Configurations of the complex
	void setConfig(std::map<std::string, std::string>&);
	void setDistanceMatrix(std::vector<std::vector<double>>* _distMatrix);
	void setIncidenceMatrix(std::vector<std::vector<bool>>* _incidenceMatrix);

	void setEnclosingRadius(double);
	static simplexBase* newSimplex(const std::string &, std::map<std::string, std::string>&);

	//Stream evaluator - this uses a function to determine if points should be inserted into the complex
	bool (*streamEval) (std::vector<double>&, std::vector<std::vector<double>>&);
    bool streamEvaluator(std::vector<double>&, std::vector<std::vector<double>>&);
	void setStreamEvaluator(bool (*f) (std::vector<double>&, std::vector<std::vector<double>>&));

	//virtual interface functions
	virtual ~simplexBase();
	virtual double getSize();
	virtual int simplexCount();
	virtual int vertexCount();
	virtual void outputComplex();

	virtual void insert();
	virtual bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&);
	virtual bool insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&);

	virtual bool find(std::vector<unsigned>);
	virtual bool find(std::set<unsigned>);

	virtual void deleteIterative(int);
	virtual void deleteIndexRecurse(int);  				// A wrapper for the actual deleteIndexRecurse method.

	virtual void prepareCofacets(int);
	virtual void prepareFacets(int);

	virtual void expandDimensions(int);
	virtual std::vector<templateNode_P> expandDimension(std::vector<templateNode_P> edges);

	//virtual void reduceComplex();



	virtual std::vector<templateNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<templateNode_P, templateNode_P>&, bool);
	virtual std::vector<nodeType*> getAllCofacets(templateNode_P, const std::unordered_map<long long, templateNode_P>&, bool);
	virtual std::vector<nodeType*> getAllCofacets(templateNode_P);
    virtual std::vector<templateNode_P> getAllCofacets(const std::set<unsigned>&);

	virtual std::vector<templateNode_P> getAllDelaunayCofacets(templateNode_P);
	virtual std::vector<nodeType*> getAllDelaunayCofacets_basePointer(templateNode_P);
	virtual std::vector<templateNode_P> getAllDelaunayCofacets(templateNode_P simp, std::unordered_map<templateNode_P,templateNode_P> pivotPairs,bool emergent);
	
	virtual std::vector<nodeType*> getAllFacets(nodeType*);
	virtual std::vector<nodeType*> getAllFacets(templateNode_P);
	virtual std::vector<templateNode_P> getAllFacets_P(templateNode_P);


	virtual std::set<templateNode_P, cmpByWeight<templateNode_P>> getDimEdges(int);
	virtual std::set<templateNode_P, cmpByWeight<templateNode_P>> getdelaunayDimEdges(int);
	virtual std::vector<std::set<templateNode_P, cmpByWeight<templateNode_P>>> getAllEdges();
	virtual std::vector<templateNode_P> expanddelaunayDimension(int);

};
