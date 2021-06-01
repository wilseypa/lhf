 
#pragma once
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>
#include "simplexArrayList.hpp" 


template <typename nodeType>
class witnessComplex : public simplexArrayList<nodeType>{
	typedef std::shared_ptr<nodeType> templateNode_P;
	
	private:

	public:
		witnessComplex(double, double);

		//virtual interface functions
		//double getSize();
		//void insert();
		//bool find(std::set<unsigned>);
		//int simplexCount();
		//int vertexCount();
		//void prepareCofacets(int);
		//void prepareFacets(int);
		//std::vector<templateNode_P> getAllCofacets(const std::set<unsigned>&, double, const std::unordered_map<templateNode_P, templateNode_P>&, bool = true);
		//bool deletion(std::set<unsigned>);
		//void expandDimensions(int);

		~witnessComplex();
		
};

template class witnessComplex<simplexNode>;
template class witnessComplex<alphaNode>;
 
