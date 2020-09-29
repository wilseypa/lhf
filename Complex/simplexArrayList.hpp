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

class simplexArrayList : public simplexBase{
	private:
		binomialTable bin;
		std::unordered_map<long long, simplexNode_P> indexConverter;

		long long simplexHash(const std::set<unsigned>&);
		unsigned maxVertex(long long, unsigned, unsigned, unsigned);
		std::vector<unsigned> getVertices(long long, int, unsigned);
	public:
		simplexArrayList(double, double);
		double findWeight(std::set<unsigned>);
		std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> recurseReduce(simplexNode_P, std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>);

		//virtual interface functions
		double getSize();
		void insert();
		bool find(std::set<unsigned>);
		int simplexCount();
		int vertexCount();
		void prepareCofacets(int);
		std::vector<simplexNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<simplexNode_P, simplexNode_P>&, bool);
		std::vector<simplexNode_P> getAllCofacets(simplexNode_P, const std::unordered_map<long long, simplexNode_P>&, bool = true);
		std::vector<simplexNode_P> getAllCofacets(simplexNode_P);
		bool deletion(std::set<unsigned>);
		void expandDimensions(int);
		void reduceComplex();
		~simplexArrayList();
};

