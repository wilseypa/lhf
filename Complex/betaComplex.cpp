#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "betaComplex.hpp"
#include <fstream>

template <typename nodeType>
betaComplex<nodeType>::betaComplex(double maxE, double maxD) : alphaComplex<nodeType>::alphaComplex(maxE, maxD)
{
	// simplexArrayList<nodeType>::simplexArrayList(0,0)

	std::cout << "Constructed Beta Complex!" << std::endl;
	this->simplexType = "betaComplex";
	this->maxEpsilon = maxE;
	this->maxDimension = maxD;
}

template <typename nodeType>
betaComplex<nodeType>::~betaComplex()
{
	this->simplexList.clear();
}

template <typename nodeType>
void betaComplex<nodeType>::buildBetaComplex(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData, double beta, std::string betaMode)
{
	unsigned maxDimension = (*dsimplexmesh.begin()).size();
	this->bin = binomialTable(npts, this->maxDimension + 1);
	for (int i = 0; i <= this->maxDimension; i++)
		this->simplexList.push_back({});

	std::ofstream out("dsimplexmesh" + std::to_string(beta) + ".csv");
	for (auto simplex : dsimplexmesh)
	{
		for (auto x : simplex)
			out << x << ",";
		out << "\n";
		unsigned int pow_set_size = pow(2, simplex.size());
		for (int counter = 1; counter < pow_set_size; counter++)
		{
			double weight = 0;
			std::set<unsigned> gensimp;
			for (int j = 0; j < simplex.size(); j++)
			{
				if (counter & (1 << j))
				{
					unsigned indnew;
					indnew = *(std::next(simplex.begin(), j));
					for (auto x : gensimp)
					{
						if (x < indnew)
						{
							if (weight < (*(this->distMatrix))[x][indnew])
								weight = (*(this->distMatrix))[x][indnew];
						}
						else if (weight < (*(this->distMatrix))[indnew][x])
							weight = (*(this->distMatrix))[indnew][x];
					}
					gensimp.insert(indnew);
				}
			}

			double weight1 = weight;
			if (gensimp.size() > 1)
				weight = utils::circumRadius(gensimp, this->distMatrix);
			else
				weight = weight / 2;
			std::shared_ptr<nodeType> tot = std::make_shared<nodeType>(nodeType(gensimp, weight1));
			if (gensimp.size() == 1)
				tot->hash = *(gensimp.begin());

			else
				tot->hash = this->simplexHash(gensimp);

			this->simplexList[gensimp.size() - 1].insert(tot);
			gensimp.clear();
		}
	}
	out.close();

	int di = 0;
	for (auto x : this->simplexList)
		std::cout << "Count of " << di++ << "-simplex ::" << x.size() << "\n";
	return;
}

template <typename nodeType>
bool betaComplex<nodeType>::checkGabriel(std::vector<double> point, std::vector<unsigned> dsimplex, std::vector<std::vector<double>> &inputData, double beta)
{
	bool intersectionCircle = false;
	bool intersectionLune = false;
	if (beta < 0)
		exit(0);
	else if (beta == 0)
		return false;
	else if (beta < 1)
		intersectionCircle = true;
	else
		intersectionLune = true;

	auto betacentersandradii = utils::calculateBetaCentersandRadius(dsimplex, inputData, this->distMatrix, beta);
	int i = 0;
	for (auto bc : betacentersandradii.first)
	{
		double distance = utils::vectors_distance(point, bc);
		if (intersectionCircle && distance > betacentersandradii.second[i])
			return false;
		if (!intersectionCircle && distance < betacentersandradii.second[i])
			return true;
		i++;
	}
	if (intersectionCircle)
		return true;
	else
		return false;
}

void buildFilteration(std::vector<std::vector<unsigned>> dsimplexmesh, int npts, std::vector<std::vector<double>> inputData, double beta, kdTree tree)
{
	std::cout << "Not implemented" << std::endl;
	return;
}

// Explicit Template Class Instantiation
template class betaComplex<simplexNode>;
template class betaComplex<alphaNode>;
template class betaComplex<witnessNode>;
