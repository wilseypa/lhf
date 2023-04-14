#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "alphaComplex.hpp"
#include <fstream>


template<typename nodeType>
alphaComplex<nodeType>::alphaComplex(double maxE, double maxD) : simplexArrayList<nodeType>::simplexArrayList(0,0) {
	/**
	    alphaComplex(double maxE, double maxD)
	 
		@brief Initializes the alpha complex
		@tparam nodeType The data type of the simplex node.
		@param maxE The max epsilon limit for complex construction.
		@param maxD The max dimension limit for complex construction.
	*/
	
	std::cout << "Constructed Alpha Complex!" << std::endl;
	
	this->simplexType = "alphaComplex";
	this->maxEpsilon = maxE;
	this->maxDimension = maxD;
}

template<typename nodeType>
alphaComplex<nodeType>::~alphaComplex(){
	/**
	    ~alphaComplex()
	 
		@brief Destructs the alpha complex and frees used memory.
		@tparam nodeType The data type of the simplex node.
	*/
	this->simplexList.clear();
}


template<typename nodeType>
std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> alphaComplex<nodeType>::getdelaunayDimEdges(int dim){
	/**
		getdelaunayDimEdges(int dim)
		 
		Maintained by Anurag
		
		@brief Get Delaunay Dimensional Edges
		@tparam nodeType The data type of the simplex node (Alpha only).
		@param dim tbd
		@return tbd
	*/
	if(dim==0)
	for(int i=0; i <= this->maxDimension; i++)
		this->simplexList.push_back({});
	if(this->simplexList[dim].size()!=0)
		return this->simplexList[dim];
#pragma omp parallel for
	for (int i = 0; i < this->dsimplexmesh.size(); ++i)
	{
		auto simplex = this->dsimplexmesh[i];
		sort(simplex.begin(), simplex.end());
		unsigned int pow_set_size = pow(2, simplex.size());
		std::set<unsigned> gensimp;
		for (int counter = 1; counter < pow_set_size; counter++)
		{
			if (__builtin_popcount(counter)!=dim+1)
				continue;
			double weight = 0;
			for (int j = 0; j < simplex.size(); j++)
			{
				if (counter & (1 << j))
				{
					unsigned indnew = simplex[j];
					for (auto x : gensimp)
					{
						if (weight < (*(this->distMatrix))[x][indnew])
							weight = (*(this->distMatrix))[x][indnew];
					}
					gensimp.insert(indnew);
				}
			}
			std::shared_ptr<nodeType> tot = std::make_shared<nodeType>(nodeType(gensimp, weight));
			if (this->simplexList[gensimp.size() - 1].find(tot) == this->simplexList[gensimp.size() - 1].end())
			{
				tot->hash = gensimp.size() > 1 ? this->simplexHash(gensimp) :  *(gensimp.begin());
#pragma omp critical
				{
					this->simplexList[gensimp.size() - 1].insert(tot);
				}
			}
			gensimp.clear();
		}
	}
	return this->simplexList[dim];
}

template <typename nodeType>
std::vector<std::shared_ptr<nodeType>> alphaComplex<nodeType>::expanddelaunayDimension(int dim) {
	/**
		expanddelaunayDimension(int dim)
		 
		Maintained by Anurag
		
		@brief Expand and get next dimension of edges
		@tparam nodeType The data type of the simplex node.
		@param dim tbd
		@return tbd
	*/
	this->simplexList[dim-1].clear();
	std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> set_simplexes=getdelaunayDimEdges(dim);
	getdelaunayDimEdges(dim+1);
	std::vector<std::shared_ptr<nodeType>> ret(set_simplexes.begin(),set_simplexes.end());
	return ret;
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> alphaComplex<nodeType>::getAllDelaunayCofacets(std::shared_ptr<nodeType> simp, std::unordered_map<std::shared_ptr<nodeType>,std::shared_ptr<nodeType>> pivotPairs,bool emergent){	
	/**
		getAllDelaunayCofacets(std::shared_ptr<nodeType> simp, std::unordered_map<std::shared_ptr<nodeType>,std::shared_ptr<nodeType>> pivotPairs, bool emergent)
	
		@brief For each pivot, get column of the pivot
		@tparam nodeType The data type of the simplex node.
		@param dim tbd
		@return tbd
	*/

	std::vector<std::shared_ptr<nodeType>> ret;
	unsigned dimension  = simp->simplex.size();
	
	for(auto iter = this->simplexList[dimension].rbegin();iter!=this->simplexList[dimension].rend();++iter){
		auto simplex = *iter;
		std::vector<unsigned> :: iterator it;
		std::vector<unsigned> v(simplex->simplex.size());
		
		it = std::set_intersection(simp->simplex.begin(),simp->simplex.end(),simplex->simplex.begin(),simplex->simplex.end(),v.begin());
		v.resize(it-v.begin());
		
		if(v.size() == simp->simplex.size()){
			ret.push_back(simplex);
			if(emergent&&simplex->weight==simp->weight){
				if(pivotPairs.find(simplex) == pivotPairs.end()) return ret;
					emergent = false;
			}
		}
	}

	return ret;
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> alphaComplex<nodeType>::getAllDelaunayCofacets(std::shared_ptr<nodeType> simp){
	/**
		getAllDelaunayCofacets(std::shared_ptr<nodeType> simp)
		
		
		@brief Get Delaunay Cofacets
		@tparam nodeType The data type of the simplex node.
		@param simp tbd
		@return tbd
	*/
	std::vector<std::shared_ptr<nodeType>> ret;
	unsigned dimension  = simp->simplex.size();
	for(auto iter = this->simplexList[dimension].rbegin();iter!=this->simplexList[dimension].rend();++iter){
		auto simplex = *iter;
		std::vector<unsigned> :: iterator it;
		std::vector<unsigned> v(simplex->simplex.size());
		
		it = std::set_intersection(simp->simplex.begin(),simp->simplex.end(),simplex->simplex.begin(),simplex->simplex.end(),v.begin());
		v.resize(it-v.begin());
		
		if(v.size() == simp->simplex.size())
			ret.push_back(simplex);
	}

	return ret;
}

template<typename nodeType>
std::vector<nodeType*> alphaComplex<nodeType>::getAllDelaunayCofacets_basePointer(std::shared_ptr<nodeType> simp){
	std::vector<nodeType*> ret;
	unsigned dimension  = simp->simplex.size();

	for(auto iter = this->simplexList[dimension].rbegin();iter!=this->simplexList[dimension].rend();++iter){
		auto simplex = *iter;
		std::vector<unsigned> :: iterator it;
		std::vector<unsigned> v(simplex->simplex.size());
		
		it = std::set_intersection(simp->simplex.begin(),simp->simplex.end(),simplex->simplex.begin(),simplex->simplex.end(),v.begin());
		v.resize(it-v.begin());


		if(v.size() == simp->simplex.size()){
			nodeType* x = new nodeType(simplex->simplex,simplex->weight);
			x->hash = simplex->hash;
			ret.push_back(x);}
	}
	return ret;
}

template<typename nodeType>

bool  alphaComplex<nodeType>::checkGabriel(std::vector<double> point, std::vector<unsigned> dsimplex ,std::vector<std::vector<double>> &inputData, double beta)
{
	bool intersectionCircle= false;
	bool intersectionLune = false;
	if(beta <0)
		exit(0);
	else if(beta==0)
		return false;
	else if(beta <1)
		intersectionCircle = true;
	else
		intersectionLune = true;

	auto betacentersandradii = utils::calculateBetaCentersandRadius(dsimplex ,inputData,this->distMatrix,beta);
	int i=0;
    	for(auto bc :betacentersandradii.first){
			double distance = utils::vectors_distance(point,bc);
			if(intersectionCircle&&distance > betacentersandradii.second[i])
				return false;
			if(!intersectionCircle&&distance < betacentersandradii.second[i])
				return true;
       i++;
		}
	if(intersectionCircle)
		return true;
	else
		return false;
}


template<>
void alphaComplex<alphaNode>::buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsmiplexmesh, int npts, std::vector<std::vector<double>> inputData){
	/**
		(alphaNode) buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)
		
		Maintained by Nick for comparing to ECC
		
		@brief Build the alpha complex from delaunay triangulation
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts 
	*/
	unsigned maxDimension = dsimplexmesh[0].size()-1;
	this->bin = binomialTable(npts, this->maxDimension+1);
	
	for(int i=0; i <= maxDimension+1; i++)
		this->simplexList.push_back({});

	for (int i = 0; i < dsimplexmesh.size(); ++i)
	{
		auto simplex = dsimplexmesh[i];
		sort(simplex.begin(), simplex.end());
		
		unsigned int pow_set_size = pow(2, simplex.size());
		std::set<unsigned> gensimp;
		
		for (int counter = 1; counter < pow_set_size; counter++)
		{
			double weight = 0;
			for (int j = 0; j < simplex.size(); j++)
			{
				if (counter & (1 << j))
				{
					unsigned indnew = simplex[j];
					for (auto x : gensimp)
					{
						if (weight < (*(this->distMatrix))[x][indnew])
							weight = (*(this->distMatrix))[x][indnew];
					}
					gensimp.insert(indnew);
				}
			}
			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp, weight));
			
			if(gensimp.size()==1)
				tot->hash = *(gensimp.begin());
			        
			else
				tot->hash = this->simplexHash(gensimp);

			this->simplexList[gensimp.size()-1].insert(tot);
			
			gensimp.clear();
		}
	}

	int di=0;
	for( auto x : this->simplexList)
		std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";
	return;
}

template<typename nodeType>
void alphaComplex<nodeType>::buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData){
	/**
		(witnessNode/simplexNode) buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)
		
		Maintained by Nick for comparing to ECC
		
		@brief Build the alpha complex from delaunay triangulation; currently a stub for witnessNode/simplexNode
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts 
	*/
    std::cout<< "alphaComplex<witnessNode, simplexNode>::buildWeightedAlphaComplex Not Implemented" << std::endl;
	return;
}

template<>
void alphaComplex<alphaNode>::buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData){
	/**
		(alphaNode) buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)
		
		Maintained by Anurag
		
		@brief Build the alpha complex from delaunay triangulation
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts 
	*/
	
	
	unsigned maxDimension = dsimplexmesh[0].size()-1;
	this->bin = binomialTable(npts, this->maxDimension+1);
    
	for(int i=0; i <= this->maxDimension; i++)
		this->simplexList.push_back({});

#pragma omp parallel for
	for (int i = 0; i < dsimplexmesh.size(); ++i)
	{
		auto simplex = dsimplexmesh[i];
		sort(simplex.begin(), simplex.end());
        
		unsigned int pow_set_size = pow(2, simplex.size());
		std::set<unsigned> gensimp;
		
		for (int counter = 1; counter < pow_set_size; counter++)
		{
			if (__builtin_popcount(counter)>this->maxDimension+1)
				continue;
			double weight = 0;
			for (int j = 0; j < simplex.size(); j++)
			{
				if (counter & (1 << j))
				{
					unsigned indnew = simplex[j];
					for (auto x : gensimp)
					{
						if (weight < (*(this->distMatrix))[x][indnew])
							weight = (*(this->distMatrix))[x][indnew];
					}
					gensimp.insert(indnew);
				}
			}
			std::shared_ptr<alphaNode> tot = std::make_shared<alphaNode>(alphaNode(gensimp, weight));
			
			if (this->simplexList[gensimp.size() - 1].find(tot) == this->simplexList[gensimp.size() - 1].end())
			{
				if (gensimp.size() > 2)
				{
					tot->circumCenter = utils::circumCenter(gensimp, inputData);
					tot->circumRadius = sqrt(utils::circumRadius(gensimp, this->distMatrix));
					tot->hash = this->simplexHash(gensimp);
				}
				else if (gensimp.size() == 2)
				{
					auto first = gensimp.begin();
					std::vector<double> R;
					std::vector<double> A = inputData[*first];
      				std::advance(first, 1);
					std::vector<double> B = inputData[*(++first)];
	   				std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
					tot->circumCenter = R;
					tot->circumRadius = sqrt(utils::circumRadius(gensimp, this->distMatrix));
					tot->hash = this->simplexHash(gensimp);
				}
				else
				{
					tot->circumRadius = weight / 2;
					tot->circumCenter = inputData[*(gensimp.begin())];
					tot->hash = *(gensimp.begin());
				}
#pragma omp critical
				{
					this->simplexList[gensimp.size() - 1].insert(tot);
				}
			}
			gensimp.clear();
		}
	}

	int di=0;
	for( auto x : this->simplexList)
		std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";
	return;
}


 
template<typename nodeType>
void alphaComplex<nodeType>::buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData){
    /**
		(witnessNode/simplexNode) buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)
		
		@brief Build the alpha complex; currently a stub for witnessNode/simplexNode
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts 
	*/
	std::cout<< "alphaComplex<witnessNode, simplexNode>::buildAlphaComplex() Not Implemented" << std::endl;
	return;
}

//Explicit Template Class Instantiation
template class alphaComplex<simplexNode>;
template class alphaComplex<alphaNode>;
template class alphaComplex<witnessNode>;

