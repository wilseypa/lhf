#include<iostream>
#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include "../Utils/utils.hpp"
#include "../Utils/readInput.hpp"
#include <bits/stdc++.h>
#include<omp.h>
using namespace std;

struct polytope{
	set<unsigned> polytopeIndices;
	set<unsigned> cofaceIndices;
	double weight;
};

class polytopalComplex{
		vector<set<polytope>> polytopalArrayList;
		vector<vector<double>> distanceMatrix;
	public:
		polytopalComplex(vector<vector<double>> &);
		std::vector<std::vector<unsigned>> qdelaunay_o(std::vector<std::vector<double>> &);
		void populateDistanceMatrix(std::vector<std::vector<double>> &);
		double getweight(vector<unsigned>);
		vector<unsigned> intersection(vector<unsigned> polytop,vector<unsigned> simplex);
		std::vector<std::vector<unsigned>> generatesimplexfacets(std::vector<unsigned> c, unsigned k);
		std::vector<std::vector<unsigned>> hullfromtriangulation(std::vector<std::vector<unsigned>> simplices);
		vector<vector<unsigned>>  generateConvexFaces(vector<vector<unsigned>> convexFaces,std::vector<std::vector<unsigned>> hull, pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<unsigned> projectionfacet,vector<double> pp);
		vector<vector<double>> projectpointsOnUnitdSphere(vector<vector<double>> &inputData,vector<double> &center);
		pair<vector<unsigned>,unsigned> findOppositeSimplex(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<vector<unsigned>> hull);
		pair<vector<vector<double>>,vector<vector<double>>> transformingConvexPolytopeForConvexDecomposition(vector<vector<double>> &inputData,vector<unsigned> projectionSimplex);
		vector<vector<double>> projectOnSimplexPlane(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<unsigned> hull,vector<double> ppoint);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> iterativeconvexization(std::vector<std::vector<unsigned>> simplices,int dim,vector<vector<double>> pointscoord);
		pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> optimalconvexization(std::vector<std::vector<unsigned>> simplices,vector<vector<double>> pointscoord,set<unsigned> originalpoints);
		pair<vector<unsigned>,set<unsigned>> getnextsimplex(std::vector<std::vector<unsigned>> simplices,pair<vector<unsigned>,set<unsigned>> cSimpAddressedpoints);
		pair<pair<vector<unsigned>,vector<vector<unsigned>>>,double> mergeneighbors(vector<unsigned> polytop,std::vector<std::vector<unsigned>> simplices,vector<vector<double>>pointscoord,vector<vector<unsigned>> delaunaypart,set<unsigned> points);
		pair<vector<vector<unsigned>>,vector<vector<unsigned>>> neighbours(vector<unsigned>  polytop,std::vector<std::vector<unsigned>> simplices,int dim1,vector<vector<unsigned>> delaunaypart);
};
