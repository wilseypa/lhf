 
#pragma once
#include "simplexBase.hpp"
#include <set>
#include <unordered_map>
#include "alphaComplex.hpp"
#include "simplexArrayList.hpp"  
#include "utils.hpp"

template <typename nodeType>
class betaComplex : public alphaComplex<nodeType> {
	typedef std::shared_ptr<nodeType> templateNode_P;
	
	private:
	      
		struct c_unique {
			unsigned current;
			c_unique() {current=0;}
			unsigned operator()() {return ++current;}
		} UniqueNumber;
		
	public:
		betaComplex(double, double);
		~betaComplex();
		
		void buildBetaComplex(std::vector<std::vector<unsigned>>dsimplexmesh,int pts, std::vector<std::vector<double>> inData,double beta,std::string betaMode);
		void buildFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData,double beta,kdTree tree);
		bool  checkGabriel(std::vector<double> , std::vector<unsigned>,std::vector<std::vector<double>>& , double );
		

};
