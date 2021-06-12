#pragma once
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>

// Header file for simplexTree class - see simplexTree.cpp for descriptions

class binomialTable{
	private:
		std::vector<std::vector<long long>> v;
	public:
		binomialTable(unsigned n, unsigned k);
		long long binom(unsigned n, unsigned k);
};

template <typename nodeType>
class simplexArrayList : public simplexBase<nodeType>{
	
	typedef std::shared_ptr<nodeType> templateNode_P;
	
	private:
		std::unordered_map<long long, templateNode_P> indexConverter;

	public:
		binomialTable bin;
		simplexArrayList(double, double);
		double findWeight(std::set<unsigned>);

		long long simplexHash(const std::set<unsigned>&);
		unsigned maxVertex(long long, unsigned, unsigned, unsigned);
		std::set<unsigned> getVertices(long long, int, unsigned);

		void initBinom();
		std::vector<nodeType*> getAllCofacets(templateNode_P, const std::unordered_map<long long, templateNode_P>&, bool = true, bool = true, unsigned = 0);
		std::vector<nodeType*> getAllCofacets(templateNode_P);
		std::vector<nodeType*> getAllFacets(nodeType*, bool = true, unsigned = 0);
		std::vector<nodeType*> getAllFacets(templateNode_P, bool = true, unsigned = 0);
		std::vector<templateNode_P> getAllFacets_P(templateNode_P);
		std::vector<templateNode_P> expandDimension(std::vector<templateNode_P>, bool = true, unsigned = 0);

		//virtual interface functions
		double getSize();
		void insert();
		bool find(std::set<unsigned>);
		int simplexCount();
		int vertexCount();
		void prepareCofacets(int);
		void prepareFacets(int);
		std::vector<templateNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<templateNode_P, templateNode_P>&, bool = true);
		bool deletion(std::set<unsigned>);
		void expandDimensions(int);

		~simplexArrayList();
};
