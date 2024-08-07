#pragma once

// Header file for betaSkeletonBasedComplexPipe class - see betaSkeletonBasedComplexPipe.cpp for descriptions
#include <map>
#include "basePipe.hpp"
#include "kdTree.hpp"

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"

using orgQhull::Coordinates;
using orgQhull::PointCoordinates;
using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullFacetSet;
using orgQhull::QhullFacetSetIterator;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullPointsIterator;
using orgQhull::QhullQh;
using orgQhull::QhullUser;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;
using orgQhull::QhullVertexSet;
using orgQhull::QhullVertexSetIterator;
using orgQhull::RboxPoints;

template <typename nodeType>
class betaSkeletonBasedComplex : public basePipe<nodeType>
{
private:
	double beta;
	std::string betaMode;
	double enclosingRadius;
	int dim;
	std::string betaMesh;
	double epsilon;
	// For generation of combinations n choose r
	struct c_unique
	{
		unsigned current;
		c_unique() { current = 0; }
		unsigned operator()() { return ++current; }
	} UniqueNumber;

public:
	betaSkeletonBasedComplex();
	void runPipe(pipePacket<nodeType> &inData);
	bool checkInsertDsimplex(std::vector<unsigned> dsimplex, pipePacket<nodeType> &inData, double beta, double averageDistance, kdTree tree);
	bool checkCC_Simplex_Inclusion(std::vector<unsigned> simplex, std::vector<std::vector<double>> inputData, std::vector<double> circumCenter);
	int getoppvertex(std::vector<unsigned> simplex, std::vector<std::vector<double>> inputData, std::vector<double> circumCenter);
	unsigned selectCenter(std::vector<double> hpcofffaces, std::vector<std::vector<double>> betaCenters, std::vector<double> otherPoint);
	std::vector<std::vector<int>> qdelaunay_o(const Qhull &qhull);
	std::vector<std::vector<int>> qconvex_o(const Qhull &qhull);
	bool configPipe(std::map<std::string, std::string> &configMap);
	void outputData(pipePacket<nodeType> &);
};
