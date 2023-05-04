#include<iostream>
#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include "../Utils/utils.hpp"
#include "../Utils/readInput.hpp"
#include <bits/stdc++.h>

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
};
