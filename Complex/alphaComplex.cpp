#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "alphaComplex.hpp"
#include "omp.h"

template <typename nodeType>
alphaComplex<nodeType>::alphaComplex(double maxE, double maxD) : simplexArrayList<nodeType>::simplexArrayList(0, 0)
{
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

template <typename nodeType>
alphaComplex<nodeType>::~alphaComplex()
{
	/**
		~alphaComplex()

		@brief Destructs the alpha complex and frees used memory.
		@tparam nodeType The data type of the simplex node.
	*/
	this->simplexList.clear();
}

template <typename nodeType>
std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> alphaComplex<nodeType>::getDimEdges(int dim)
{
	/**
		getDimEdges(int dim)

		Maintained by Anurag

		@brief Get Delaunay Dimensional Edges
		@tparam nodeType The data type of the simplex node (Alpha only).
		@param dim tbd
		@return tbd
	*/
	if (dim == 0)
		for (int i = 0; i <= this->maxDimension; i++)
			this->simplexList.push_back({});
	if (this->simplexList[dim].size() != 0)
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
			if (__builtin_popcount(counter) != dim + 1)
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
				tot->hash = gensimp.size() > 1 ? this->simplexHash(gensimp) : *(gensimp.begin());
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
std::vector<std::shared_ptr<nodeType>> alphaComplex<nodeType>::expanddelaunayDimension(int dim)
{
	/**
		expanddelaunayDimension(int dim)

		Maintained by Anurag

		@brief Expand and get next dimension of edges
		@tparam nodeType The data type of the simplex node.
		@param dim tbd
		@return tbd
	*/
	this->simplexList[dim - 1].clear();
	std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> set_simplexes = getDimEdges(dim);
	getDimEdges(dim + 1);
	std::vector<std::shared_ptr<nodeType>> ret(set_simplexes.begin(), set_simplexes.end());
	return ret;
}

template <typename nodeType>
std::vector<std::shared_ptr<nodeType>> alphaComplex<nodeType>::getAllCofacets(std::shared_ptr<nodeType> simp)
{
	/**
		getAllDelaunayCofacets(std::shared_ptr<nodeType> simp)


		@brief Get Delaunay Cofacets
		@tparam nodeType The data type of the simplex node.
		@param simp tbd
		@return tbd
	*/
	std::vector<std::shared_ptr<nodeType>> ret;
	unsigned dimension = simp->simplex.size();
	for (auto &iter : this->simplexList[dimension])
	{
		if (std::includes(iter->simplex.begin(), iter->simplex.end(), simp->simplex.begin(), simp->simplex.end()))
			ret.push_back(iter);
	}

	return ret;
}

template <typename nodeType>
std::vector<nodeType *> alphaComplex<nodeType>::getAllCofacets_basePointer(std::shared_ptr<nodeType> simp)
{
	std::vector<nodeType *> ret;
	unsigned dimension = simp->simplex.size();
	for (auto &iter : this->simplexList[dimension])
	{
		if (std::includes(iter->simplex.begin(), iter->simplex.end(), simp->simplex.begin(), simp->simplex.end()))
		{
			nodeType *x = new nodeType(iter->simplex, iter->weight);
			x->hash = iter->hash;
			ret.push_back(x);
		}
	}
	return ret;
}

template <>
void alphaComplex<alphaNode>::buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsmiplexmesh, int npts, std::vector<std::vector<double>> inputData)
{
	/**
		(alphaNode) buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)

		Maintained by Nick for comparing to ECC

		@brief Build the alpha complex from delaunay triangulation
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts
	*/
	unsigned maxDimension = dsimplexmesh[0].size() - 1;
	this->bin = binomialTable(npts, this->maxDimension + 1);

	for (int i = 0; i <= this->maxDimension; i++)
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
			if (__builtin_popcount(counter) > this->maxDimension + 1)
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
				if (gensimp.size() == 1)
					tot->hash = *(gensimp.begin());

				else
					tot->hash = this->simplexHash(gensimp);
#pragma omp critical
				{
					this->simplexList[gensimp.size() - 1].insert(tot);
				}
			}

			gensimp.clear();
		}
	}

	int di = 0;
	for (auto x : this->simplexList)
		std::cout << "Count of " << di++ << "-simplex ::" << x.size() << "\n";
	return;
}

template <typename nodeType>
void alphaComplex<nodeType>::buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)
{
	/**
		(witnessNode/simplexNode) buildWeightedAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)

		Maintained by Nick for comparing to ECC

		@brief Build the alpha complex from delaunay triangulation; currently a stub for witnessNode/simplexNode
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts
	*/
	std::cout << "alphaComplex<witnessNode, simplexNode>::buildWeightedAlphaComplex Not Implemented" << std::endl;
	return;
}

template <>
void alphaComplex<alphaNode>::buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)
{
	/**
		(alphaNode) buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)

		Maintained by Anurag

		@brief Build the alpha complex from delaunay triangulation
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts
	*/

	unsigned maxDimension = dsimplexmesh[0].size() - 1;
	this->bin = binomialTable(npts, this->maxDimension + 1);

	for (int i = 0; i <= this->maxDimension; i++)
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
			if (__builtin_popcount(counter) > this->maxDimension + 1)
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
					std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R), [](double e1, double e2)
								   { return ((e1 + e2) / 2); });
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

	int di = 0;
	for (auto x : this->simplexList)
		std::cout << "Count of " << di++ << "-simplex ::" << x.size() << "\n";

	neighbourhood = std::vector<std::vector<int>>(npts, std::vector<int>(npts, 0));
	for (const auto &x : this->simplexList[1]) // Use edges to build neighbourhood relationship
	{
		neighbourhood[*x->simplex.begin()][*(--x->simplex.end())] = 1;
		neighbourhood[*(--x->simplex.end())][*x->simplex.begin()] = 1;
	}
	return;
}

template <typename nodeType>
void alphaComplex<nodeType>::buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)
{
	/**
		(witnessNode/simplexNode) buildAlphaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData)

		@brief Build the alpha complex; currently a stub for witnessNode/simplexNode
		@tparam nodeType The data type of the simplex node.
		@param dsimplexmesh the set of d-triangles for the simplex mesh
		@param npts
	*/
	std::cout << "alphaComplex<witnessNode, simplexNode>::buildAlphaComplex() Not Implemented" << std::endl;
	return;
}

// Explicit Template Class Instantiation
template class alphaComplex<simplexNode>;
template class alphaComplex<alphaNode>;
template class alphaComplex<witnessNode>;
