/*
 * betaSubSkeletonComplexPipe hpp + cpp extend the basePipe class for calculating the
 * beta Skeleton Based Complex generation for data input
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <set>
#include <algorithm>
#include "betaSubSkeletonComplex.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include "readInput.hpp"

// basePipe constructor
template <typename nodeType>
betaSubSkeletonComplex<nodeType>::betaSubSkeletonComplex()
{
	this->pipeType = "betaSubSkeletonComplex";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void betaSubSkeletonComplex<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	// Generate Beta Skeleton Based Complex
	// Temporarily commenting this out - need to check inData.complex type
	//		If type is the graph-based simplexArrayList (inherited) then
	//			cast to gAL and run non-virtual function:

	//	((graphArrayList*)inData.complex)->graphInducedComplex(dim,inData.inputData,beta);
	std::vector<std::vector<unsigned>> dsimplexmesh;
	kdTree tree(inData.inputData, inData.inputData.size()); // KDTree for efficient nearest neighbor search
	// int dim = inData.inputData[0].size();
	int dim = this->dim;
	double distanceSum = 0;
	int count = 0;
	// Computing the average point cloud distance
	for (int x = 0; x < inData.inputData.size(); x++)
		for (int y = x + 1; y < inData.inputData.size(); y++)
		{
			distanceSum += (*((alphaComplex<nodeType> *)inData.complex)->distMatrix)[x][y];
			count = count + 1;
		}
	double averageDistance = distanceSum / count;
	count = 0;
	std::vector<std::pair<int, int>> neighborsepsilon;

	// Enumeration within epsilon ball for beta <1.
	if (this->betaMesh == "null.csv")
	{
		for (unsigned index = 0; index < inData.inputData.size(); index++)
		{
			std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); // All neighbors in epsilon-ball
			int n = neighbors.size();
			neighborsepsilon.push_back(std::make_pair(n, index));
		}
		sort(neighborsepsilon.begin(), neighborsepsilon.end());
		std::vector<int> toremove;
		for (unsigned indexold = 0; indexold < inData.inputData.size(); indexold++)
		{
			toremove.push_back(neighborsepsilon[indexold].second);
			unsigned index = neighborsepsilon[indexold].second;
			std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); // All neighbors in epsilon-ball
			int n = neighbors.size();
			std::sort(toremove.begin(), toremove.end());
			std::sort(neighbors.begin(), neighbors.end());
			std::vector<int> difference;
			std::set_difference(neighbors.begin(), neighbors.end(), toremove.begin(), toremove.end(), std::back_inserter(difference));
			n = difference.size();
			if (n >= dim)
			{
				std::vector<unsigned> dsimplex(dim);
				std::vector<unsigned> dsimplexIndexed;
				std::vector<unsigned>::iterator first = dsimplex.begin(), last = dsimplex.end();
				std::generate(first, last, UniqueNumber);
				dsimplexIndexed.push_back(index);
				for (int i = 0; i < dim; i++)
					dsimplexIndexed.push_back(difference[dsimplex[i] - 1]);
				// for each enumerated simplex evaluate for beta selecton criterea
				if (checkInsertSubDsimplex(dsimplexIndexed, inData, this->beta, averageDistance, tree))
				{
					std::sort(dsimplexIndexed.begin(), dsimplexIndexed.end());
					dsimplexmesh.push_back(dsimplexIndexed);
					count++;
				}
				while ((*first) != n - dim + 1)
				{
					std::vector<unsigned>::iterator mt = last;
					while (*(--mt) == n - (last - mt) + 1)
						;
					(*mt)++;
					while (++mt != last)
						*mt = *(mt - 1) + 1;
					std::vector<unsigned> dsimplexIndexed1;
					dsimplexIndexed1.push_back(index);
					for (int i = 0; i < dim; i++)
					{
						dsimplexIndexed1.push_back(difference[dsimplex[i] - 1]);
						//   std::cout<<difference[dsimplex[i]-1]<<" ";
					}
					// for each enumerated simplex evaluate for beta selecton criterea
					if (checkInsertSubDsimplex(dsimplexIndexed1, inData, this->beta, averageDistance, tree))
					{
						std::sort(dsimplexIndexed1.begin(), dsimplexIndexed1.end());
						dsimplexmesh.push_back(dsimplexIndexed1);
						count++;
					}
				}
			}
		}
		std::ofstream out("latestbetamesh.csv");
		for (auto x : dsimplexmesh)
		{
			int i = 0;
			for (auto y : x)
			{
				i++;
				out << y;
				if (i != x.size())
					out << ",";
			}
			out << '\n';
		}
	}
	else if (betaMesh != "null.csv")
	{
		auto rs = readInput();
		std::vector<std::vector<double>> betaMesh1 = rs.readCSV(this->betaMesh);
		std::vector<std::vector<unsigned>> betaMesh;
		for (auto x : betaMesh1)
		{
			std::vector<unsigned> temp;
			for (auto y : x)
				temp.push_back(y);
			betaMesh.push_back(temp);
		}
		for (auto x : betaMesh)
			// for each enumerated simplex evaluate for beta selecton criterea
			if (checkInsertSubDsimplex(x, inData, this->beta, averageDistance, tree))
			{
				std::sort(x.begin(), x.end());
				dsimplexmesh.push_back(x);
			}
		std::ofstream out("latestbetamesh.csv");
		for (auto x : dsimplexmesh)
		{
			int i = 0;
			for (auto y : x)
			{
				i++;
				out << y;
				if (i != x.size())
					out << ",";
			}
			out << '\n';
		}
	}

	//	((alphaComplex<nodeType>*)inData.complex)->buildBetaComplex(dsimplexmesh,inData.inputData.size(),inData.inputData,this->beta,this->betaMode);
	//	((alphaComplex<nodeType>*)inData.complex)->buildAlphaComplex(dsimplexmesh,inData.inputData.size(),inData.inputData);
	//	((alphaComplex<nodeType>*)inData.complex)->buildFilteration(dsimplexmesh,inData.inputData.size(),inData.inputData,this->beta);
	//        ((alphaComplex<nodeType>*)inData.complex)->buildBetaComplexFilteration(dsimplexmesh, inData.inputData.size(),inData.inputData, tree);
	std::vector<std::vector<bool>> incidenceMatrix(inData.inputData.size(), std::vector<bool>(inData.inputData.size(), 0));
	int countd1 = 0;
	for (auto x : dsimplexmesh)
	{
		for (int i = 0; i < x.size(); i++)
			for (int j = i + 1; j < x.size(); j++)
			{
				int origin = x[i];
				int destination = x[j];
				if (incidenceMatrix[origin][destination] != true)
					countd1++;
				incidenceMatrix[origin][destination] = true;
			}
	}

	// std::cout<<"Edges ="<<countd1<<" ";
	inData.incidenceMatrix = incidenceMatrix;
	std::ofstream file("PHdSphereDimensionWiseMeshSize.txt", std::ios_base::app);
	file << this->betaMode << "," << inData.inputData.size() << "," << inData.inputData[0].size() << "," << this->beta << "," << dsimplexmesh.size() << std::endl;
	file.close();
	this->ut.writeDebug("betaSubSkeletonComplex Pipe", "\tbetaSubSkeletonComplex Size: ");
	return;
}

template <>
bool betaSubSkeletonComplex<alphaNode>::checkInsertSubDsimplex(std::vector<unsigned> dsimplex, pipePacket<alphaNode> &inData, double beta, double averageDistance, kdTree tree)
{
	double maxEdge = 0;
	for (auto x : dsimplex)
		for (auto y : dsimplex)
			if (maxEdge < inData.distMatrix[x][y])
				maxEdge = inData.distMatrix[x][y];

	if (maxEdge > this->epsilon)
		return false;

	std::vector<size_t> neighborsfinalLune;
	std::vector<size_t> neighborsfinalCircle;

	bool intersectionCircle = false;
	bool intersectionLune = false;
	bool unionCircle = false;

	if (beta < 0)
		exit(0);
	else if (beta == 0)
		return true;
	else if (beta < 1)
		intersectionCircle = true;
	else if (beta == 1)
	{
		intersectionLune = true;
		intersectionCircle = true;
		unionCircle = true;
	}
	else if (beta > 1)
	{
		if (this->betaMode == "betaHighLune")
			intersectionLune = true;
		else
			unionCircle = true;
	}
	std::vector<double> circumCenter;
	std::set<unsigned> simplex(dsimplex.begin(), dsimplex.end());

	if (simplex.size() > 2)
		circumCenter = utils::circumCenter(simplex, inData.inputData);
	else if (simplex.size() == 2)
	{
		auto first = simplex.begin();
		std::vector<double> R;
		std::vector<double> A = inData.inputData[*first];
		std::advance(first, 1);
		std::vector<double> B = inData.inputData[*first];
		std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R), [](double e1, double e2)
					   { return ((e1 + e2) / 2); });
		circumCenter = R;
	}
	double circumRadius;
	circumRadius = pow(utils::vectors_distance(circumCenter, inData.inputData[dsimplex[0]]), 2);

	if (this->betaMode == "betaHighCircle")
	{
		if (beta < 1)
			beta = 1 / beta;

		std::vector<size_t> neighbors;
		std::vector<std::vector<size_t>> neighborsCircleIntersection;
		std::vector<std::vector<double>> mat;

		for (auto x : simplex)
			mat.push_back(inData.inputData[x]);
		mat.push_back(circumCenter);
		auto hpcoff = utils::nullSpaceOfMatrix(simplex, inData.inputData, circumCenter, sqrt(circumRadius), true);

		std::vector<std::vector<double>> refbetaCenters;
		refbetaCenters = utils::betaCentersCalculation(hpcoff.first, beta, sqrt(circumRadius), circumCenter);

		refbetaCenters = utils::computePCAInverse(mat, refbetaCenters, hpcoff.second);
		double betaRadius = utils::vectors_distance(refbetaCenters[0], inData.inputData[dsimplex[0]]);
		std::vector<size_t> neighbors1 = tree.neighborhoodIndices(refbetaCenters[0], betaRadius); // All neighbors in epsilon-ball
		for (auto t : dsimplex)
			neighbors1.erase(std::remove(neighbors1.begin(), neighbors1.end(), t), neighbors1.end());
		std::vector<size_t> neighbors2 = tree.neighborhoodIndices(refbetaCenters[1], betaRadius); // All neighbors in epsilon-ball
		std::sort(neighbors1.begin(), neighbors1.end());
		std::sort(neighbors2.begin(), neighbors2.end());
		if (intersectionCircle == true)
		{
			std::vector<size_t> v1(std::min(neighbors1.size(), neighbors2.size()));
			std::vector<size_t>::iterator it1;
			it1 = std::set_intersection(neighbors1.begin(), neighbors1.end(), neighbors2.begin(), neighbors2.end(), v1.begin());
			v1.resize(it1 - v1.begin());
			neighborsfinalCircle = v1;
		}
		else if (unionCircle == true)
		{
			std::vector<size_t> v1(neighbors1.size() + neighbors2.size());
			std::vector<size_t>::iterator it1;
			it1 = std::set_union(neighbors1.begin(), neighbors1.end(), neighbors2.begin(), neighbors2.end(), v1.begin());
			v1.resize(it1 - v1.begin());
			neighborsfinalCircle = v1;
		}
	}

	else if (this->betaMode == "betaHighLune")
	{
		std::vector<double> circumCenterfaces;
		std::vector<double> circumCenterfaces1;

		bool first = true;
		for (auto x : simplex)
		{
			std::vector<double> betaCenter;
			for (unsigned y = 0; y < inData.inputData[0].size(); y++)
				betaCenter.push_back(beta * circumCenter[y] + (1 - beta) * inData.inputData[x][y]);
			double betaRadius = beta * sqrt(circumRadius);

			std::vector<size_t> neighbors1faces1 = tree.neighborhoodIndices(betaCenter, betaRadius); // All neighbors in epsilon-ball
			neighbors1faces1.erase(std::remove(neighbors1faces1.begin(), neighbors1faces1.end(), x), neighbors1faces1.end());
			if (!first)
			{
				std::sort(neighborsfinalLune.begin(), neighborsfinalLune.end());
				std::sort(neighbors1faces1.begin(), neighbors1faces1.end());
				if (intersectionLune == true)
				{
					std::vector<size_t> v1(std::min(neighborsfinalLune.size(), neighbors1faces1.size()));
					std::vector<size_t>::iterator it1;
					it1 = std::set_intersection(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1 - v1.begin());
					neighborsfinalLune = v1;
				}
				else
				{
					std::vector<size_t> v1(neighborsfinalLune.size() + neighbors1faces1.size());
					std::vector<size_t>::iterator it1;
					it1 = std::set_union(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1 - v1.begin());
					neighborsfinalLune = v1;
				}
			}
			else
				neighborsfinalLune = neighbors1faces1;
			first = false;
		}
	}
	if (this->betaMode == "betaHighLune" && neighborsfinalLune.size() == 0)
		return true;
	if (this->betaMode == "betaHighCircle" && neighborsfinalCircle.size() == 0)
		return true;
	return false;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool betaSubSkeletonComplex<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
{
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if (pipe != configMap.end())
	{
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if (pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();

	pipe = configMap.find("beta");
	if (pipe != configMap.end())
		this->beta = std::atof(configMap["beta"].c_str());

	pipe = configMap.find("betaMode");
	if (pipe != configMap.end())
		this->betaMode = configMap["betaMode"].c_str();

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->epsilon = std::atof(configMap["epsilon"].c_str());

	this->ut = utils(strDebug, this->outputFile);
	pipe = configMap.find("dimensions");
	if (pipe != configMap.end())
	{
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	pipe = configMap.find("betaMesh");
	if (pipe != configMap.end())
		this->betaMesh = configMap["betaMesh"].c_str();

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else
		return false;

	this->configured = true;
	this->ut.writeDebug("betaSubSkeletonComplex Pipe ", "Configured with parameters { eps: " + configMap["epsilon"] + configMap["beta"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void betaSubSkeletonComplex<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	// Output related to betaSubSkeletonComplex
	return;
}
template <typename nodeType>
bool betaSubSkeletonComplex<nodeType>::checkInsertSubDsimplex(std::vector<unsigned> dsimplex, pipePacket<nodeType> &inData, double beta, double averageDistance, kdTree tree)
{
	std::cout << "No function Defined";
	return false;
}

template class betaSubSkeletonComplex<simplexNode>;
template class betaSubSkeletonComplex<alphaNode>;
template class betaSubSkeletonComplex<witnessNode>;
