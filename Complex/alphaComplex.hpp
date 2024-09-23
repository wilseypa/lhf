
#pragma once
#include <set>
#include <unordered_map>
#include "simplexArrayList.hpp"
#include "utils.hpp"
#include "kdTree.hpp"

template <typename nodeType>
class alphaComplex : public simplexArrayList<nodeType>
{
	typedef std::shared_ptr<nodeType> templateNode_P;

public:
	alphaComplex(double, double);

	// virtual interface functions
	// double getSize();
	// void insert();
	// bool find(std::set<unsigned>);
	// int simplexCount();
	// int vertexCount();
	// void prepareCofacets(int);
	// void prepareFacets(int);
	// std::vector<templateNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<templateNode_P, templateNode_P>&, bool = true);
	// bool deletion(std::set<unsigned>);
	// void expandDimensions(int);
	~alphaComplex();

	std::vector<templateNode_P> getAllDelaunayCofacets(templateNode_P);
	std::vector<templateNode_P> getAllDelaunayCofacets(templateNode_P simp, std::unordered_map<templateNode_P, templateNode_P> pivotPairs, bool emergent);
	std::vector<nodeType *> getAllDelaunayCofacets_basePointer(templateNode_P);
	void buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int pts, std::vector<std::vector<double>> inputData);
	void buildFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData, double beta, kdTree tree);
	bool checkGabriel(std::vector<double>, std::vector<unsigned>, std::vector<std::vector<double>> &, double);
	std::vector<templateNode_P> expanddelaunayDimension(int);
	std::set<templateNode_P, cmpByWeight<templateNode_P>> getDimEdges(int);

	// Nick

	void buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsmiplexmesh, int npts, std::vector<std::vector<double>> inputData);
};
