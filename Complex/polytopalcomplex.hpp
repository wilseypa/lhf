#include <iostream>
#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include "../Utils/utils.hpp"
#include "../Utils/readInput.hpp"
#include <bits/stdc++.h>
#include <omp.h>
#include <unordered_map>
#include <cddlib/setoper.h>
#include <cddlib/cdd.h>

using namespace std;

struct polytope{
	set<unsigned> polytopeIndices;
	vector<vector<unsigned>> Deltriangulation;
	vector<vector<double>> coordinates;
	set<unsigned> cofaceIndices;
	double weight;
	bool operator==(const polytope& p) const
    {
        return polytopeIndices == p.polytopeIndices && weight == p.weight;
    }
    bool operator > (const  polytope&p) const {
		if(weight==p.weight)
			return polytopeIndices <p.polytopeIndices;
		
        return weight>p.weight;
    }
    bool operator <(const polytope&p) const {
		if(weight==p.weight)
			return polytopeIndices >p.polytopeIndices;
		
        return weight<p.weight;
    }
};
struct cmp {
    bool operator() (struct  polytope p1, struct  polytope p2) const {
		if(p1.weight==p2.weight)
			return p1.polytopeIndices >p2.polytopeIndices;
		
        return p1.weight>p2.weight;
    }
};
class MyHashFunction {
public:
 
    std::size_t operator()(const polytope& vec1) const {
	std::size_t seed = vec1.polytopeIndices.size();
	std::vector<uint32_t> vec(vec1.polytopeIndices.begin(),vec1.polytopeIndices.end());
	for(auto x : vec) {
		x = ((x >> 16) ^ x) * 0x45d9f3b;
		x = ((x >> 16) ^ x) * 0x45d9f3b;
		x = (x >> 16) ^ x;
		seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}
  return seed;
}
};

class polytopalComplex{
		vector<set<polytope,cmp>> polytopalArrayList;
		vector<vector<double>> distanceMatrix;
		vector<vector<double>> bettiTable;
		int dataSize;
		int dim;
	public:
		vector<vector<double>> getbettiTable(){
			return bettiTable;
		}
		polytopalComplex(vector<vector<double>> &);
		std::vector<std::vector<unsigned>> qdelaunay_o(std::vector<std::vector<double>> &);
		void populateDistanceMatrix(std::vector<std::vector<double>> &);
		double getweight(vector<unsigned>);
		double  distanceindex(unsigned,unsigned);
		vector<unsigned> intersection(vector<unsigned> polytop,vector<unsigned> simplex);
		std::vector<std::vector<unsigned>> generatesimplexfacets(std::vector<unsigned> c, unsigned k);
		std::vector<std::vector<unsigned>> hullfromtriangulation(std::vector<std::vector<unsigned>> &simplices);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> generateConvexFaces(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexFaces,std::vector<std::vector<unsigned>> hull, pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<unsigned> projectionfacet,vector<double> pp);
		vector<vector<double>> projectpointsOnUnitdSphere(vector<vector<double>> &inputData,vector<double> &center);
		pair<vector<unsigned>,unsigned> findOppositeSimplex(pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<vector<unsigned>> hull);
		pair<vector<vector<double>>,vector<vector<double>>> transformingConvexPolytopeForConvexDecomposition(vector<vector<double>> &inputData,vector<unsigned> projectionSimplex);
		vector<vector<double>> projectOnSimplexPlane(pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<double> ppoint,double scaledown);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> iterativeconvexization(std::vector<std::vector<unsigned>> &simplices,int dim,vector<vector<double>> &pointscoord);
		pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> optimalconvexization(std::vector<std::vector<unsigned>> &simplices,vector<vector<double>> &pointscoord,set<unsigned> &originalpoints);
		pair<vector<unsigned>,set<unsigned>> getnextsimplex(std::vector<std::vector<unsigned>> &simplices,pair<vector<unsigned>,set<unsigned>> cSimpAddressedpoints,std::vector<std::vector<unsigned>> addressedsimplices);
		pair<pair<vector<unsigned>,vector<vector<unsigned>>>,double> mergeneighbors(vector<unsigned> polytop,std::vector<std::vector<unsigned>> &simplices,vector<vector<double>>&pointscoord,vector<vector<unsigned>> delaunaypart,set<unsigned> &points);
		pair<pair<vector<unsigned>,vector<vector<unsigned>>>,double> mergeneighborsPrev	(vector<unsigned> polytop,std::vector<std::vector<unsigned>> &simplices,vector<vector<double>>&pointscoord,vector<vector<unsigned>> delaunaypart,set<unsigned> &points);
		pair<vector<vector<unsigned>>,vector<vector<unsigned>>> neighbours(vector<unsigned>  polytop,std::vector<std::vector<unsigned>> &simplices,int dim1,vector<vector<unsigned>> delaunaypart);
		pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> sortofSizeandWeight(pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> &df);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  informedConvexization(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  &convexFaces, std::vector<std::vector<unsigned>> &hull,pair<vector<vector<double>>,vector<double>> &projectionData);
		pair<vector<vector<double>>,vector<double>> projectonCenteroidSphere(vector<vector<double>>&inputData);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> sortFaces(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexfaces);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> pruneMaximalParts(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexFaces,pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexPolytopes);
		bool checkCC_Simplex_Inclusion(std::vector<unsigned> simplex,std::vector<std::vector<double> > &inputData,	std::vector<double> circumCenter);
		std::vector<std::vector<unsigned>>  transformhull(std::vector<std::vector<unsigned>> originalhull,set<unsigned> polytopeIndices); // Global to Local Indices
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  transformCF(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexFaces,set<unsigned> polytopeIndices); // local to global Indices
		double getmaxWeight(vector<vector<unsigned>> &delparts);
		void persistence();
		std::vector<polytope> persistenceByDim(std::vector<polytope> pivots,int d);
		std::vector<polytope> getCofacet(set<unsigned> cofaceIndices,int codim);
		set<polytope> generateFaces(int d);
		std::vector<polytope> getCofacetModified(set<unsigned> cofaceIndices,int codim);
		vector<vector<double>> inverseStereoGraphicProjection(vector<vector<double>> &);
		vector<vector<double>> rescaledataandcenteraroundorigin(vector<vector<double>> &data);
		vector<vector<double>> HyperplaneToVertexRepresentation(std::vector<std::vector<double>>& A,std::vector<double> &B);
		template <typename T>
		dd_MatrixPtr dd_PolyFile2Matrix_2(std::vector<std::vector<T>> A, std::vector<T> B, dd_ErrorType *Error);
		
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> approxDecomposition(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> Faces,vector<vector<double>> coordinates);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> mergeApprox(pair<vector<unsigned>,vector<vector<unsigned>>> x,pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> Neighbors,double approxscale,vector<vector<double>>);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> updateUnprocessed(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> unprocessed, pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> processed);
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> postprocess(std::vector<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>> toprocess);

};
