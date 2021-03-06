 
#pragma once
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>
#include "simplexArrayList.hpp" 
#include "utils.hpp"

template <typename nodeType>
class alphaComplex : public simplexArrayList<nodeType>{
	typedef std::shared_ptr<nodeType> templateNode_P;
	
	private:
	public:
		alphaComplex(double, double);

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
		~alphaComplex();
		
		std::vector<templateNode_P> getAllDelaunayCofacets(templateNode_P);
		void buildAlphaComplex(std::vector<std::vector<int>> dsimplexmesh, int pts,std::vector<std::vector<double>> inputData);
		void buildBetaComplex(std::set<std::vector<unsigned>>dsimplexmesh,int pts, std::vector<std::vector<double>> inData);

};
