#include<iostream>
#include "polytopalcomplex.hpp"

using namespace std;

double polytopalComplex :: getweight(vector<unsigned> simplex){
	double weight =0;
	for(int i=0;i<simplex.size();i++){
		for(int j=i+1;j<simplex.size();j++){
			if(weight<distanceMatrix[simplex[i]][simplex[j]])
				weight = distanceMatrix[simplex[i]][simplex[j]];
		}
	}
	return weight;
}

std::vector<std::vector<unsigned>>  polytopalComplex :: qdelaunay_o(vector<vector<double>> &inputData){
  int dim = inputData[0].size();
  typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>> T;
  T dt(dim);
  std::vector<T::Point> points;
  unsigned index = 0;
  std::map<T::Point, unsigned> index_of_vertex;
  for (auto i : inputData)
  {
    T::Point p(i.begin(), i.end());
    points.push_back(p);
    index_of_vertex[p] = index++;
  }
  dt.insert(points.begin(), points.end());
  points.clear();
  long unsigned int number_of_finite_full_cell = dt.number_of_finite_full_cells();
  std::vector<std::vector<unsigned>> dsimplexes(number_of_finite_full_cell, std::vector<unsigned>(dim + 1));
#pragma omp parallel for
  for (int j = 0; j <= dim; ++j)
  {
    auto i = dt.finite_full_cells_begin();
    for (long unsigned int k = 0; k < number_of_finite_full_cell; k++)
      dsimplexes[k][j] = index_of_vertex[(i++)->vertex(j)->point()];
  }
  return dsimplexes;
}

void  polytopalComplex :: populateDistanceMatrix(vector<vector<double>> &inputData){
	distanceMatrix.resize(inputData.size());
	for(auto x=0;x<inputData.size();x++){
		for(auto y=x+1;y<inputData.size();y++){
			distanceMatrix[x].push_back(utils::vectors_distance(inputData[x],inputData[y]));
		}
	}
}
std::vector<std::vector<unsigned>> generatesimplexfacets(std::vector<unsigned> c, unsigned k){
	std::vector<std::vector<unsigned>> combinations;
	unsigned long n = c.size();
    unsigned long combo = (1 << k) - 1;    
    while (combo < 1<<n) {
		std::vector<unsigned> temp;
        for (unsigned long i = 0; i < c.size(); ++i) {
			if ((combo >> i) & 1){
				temp.push_back(c[i]);
			}
		}
		combinations.push_back(temp);
        unsigned long  x = combo & -combo;
        unsigned long  y = combo + x;
        unsigned long  z = (combo & ~y);
        combo = z / x;
        combo >>= 1;
        combo |= y;
    }
    return combinations;
}
template <typename T>
bool contains(std::vector<T> first, std::vector<T> second)
{
    std::sort(first.begin(), first.end());
    std::sort(second.begin(), second.end());
    return std::includes(first.begin(), first.end(), second.begin(), second.end());
} 
std::vector<std::vector<unsigned>> hullfromtriangulation(std::vector<std::vector<unsigned>> simplices){
	std::vector<std::vector<unsigned>> hull;
	for (auto x : simplices){
		auto facets = generatesimplexfacets(x,x.size()-1);
		for (auto f : facets){
			bool valid = true;
			for (auto y : simplices){
				if(valid==false)
					break;
				if(!(x==y)){
					if(contains(y,f)){
						valid=false;
						break;
					}
				}
			}
			if(valid==true){
				std::sort(f.begin(), f.end());
				bool present = false;
				for(auto k :hull)
					if(k==f){
						present = true;
						break;
					}
				if(!present)
					hull.push_back(f);
				}
		}
	}
	return hull;
}


vector<vector<double>> projectpointsOnUnitdSphere(vector<vector<double>> &inputData,vector<double> &center){
	vector<vector<double>> updatedCoordinates;
	for(auto x : inputData){
		vector<double> point;
		double d =utils::vectors_distance(center,x);
		for(int i=0;i<inputData[0].size();i++){
			point.push_back(((((x[i]-center[i]))/d))+center[i]);
		}
		updatedCoordinates.push_back(point);
	}
	return updatedCoordinates;
}


vector<unsigned> findOppositeSimplex(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<vector<unsigned>> hull){
	vector<double> point = projectionData.second[0];
	double minweight = 9999;
	unsigned projectionSimplex;
	int i=0;
	for(auto x : hull){
		vector<double> centroid(projectionData.first[0].size(),0.0);
		for(auto y :x){
			for(int j=0;j<projectionData.first[0].size();j++){
				centroid[j] +=((projectionData.first[y][j])/projectionData.first[0].size());
			}
		}
		auto d = utils::vectors_distance(centroid,point);
		if(d<minweight){
			minweight = d;
			projectionSimplex = i; 
		}
		i++;
	}
	return hull[projectionSimplex];
}

pair<vector<vector<double>>,vector<vector<double>>> transformingConvexPolytopeForConvexDecomposition(vector<vector<double>> &inputData,vector<unsigned> projectionSimplex){
	std::vector<std::vector<double>> pts;
	std::vector<std::vector<double>> pp;
	vector<double> center(inputData[0].size(),0.0);
	bool found = false;
	for( auto b : projectionSimplex){
		pts.push_back(inputData[b]);
		for(int i=0;i<inputData[0].size();i++){
			center[i] +=((inputData[b][i])/inputData[0].size());
		}
	}
	vector<double> pp1(inputData[0].size()-1,0.0);
	vector<double> pp2(inputData[0].size()-1,0.0);
	pp1.push_back(1.0);
	pp2.push_back(-1.0);
	auto pca = utils:: computePCA(pts,pts[0].size());
	auto modifiedPoints = pca.first;
	
	modifiedPoints.push_back(pp1);
	modifiedPoints.push_back(pp2);
	auto newpoints  = utils::computePCAInverse(pts,modifiedPoints,pca.second);
	for(auto f:newpoints){
		for(auto d:f)
			cout<<d<<" ";
		cout<<"\n";
	}
	
	pp1 = newpoints[newpoints.size()-2];
	pp2 = newpoints[newpoints.size()-1];
	auto trasformedData = utils::projectHSphere(inputData,center);
	pp.push_back(pp1);
	pp.push_back(pp2);
	return make_pair(trasformedData,pp);
}
projectOnSimplexPlane(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<vector<unsigned>> hull){
	
	
	
	return ProjectedData;
}

polytopalComplex :: polytopalComplex(vector<vector<double>> &inputData){
	populateDistanceMatrix(inputData);
	auto simplices = qdelaunay_o(inputData);
	auto hull = hullfromtriangulation(simplices);
	auto weight1 = getweight(hull[0]);
	auto projectionData1 = transformingConvexPolytopeForConvexDecomposition(inputData,hull[0]);
	auto projectedonPlane1 = projectOnSimplexPlane(projectionData1,hull[0]);
	auto oppositeSimplex = findOppositeSimplex(projectionData1,hull);
	cout<<hull[0][0]<<" "<<hull[0][1]<<" "<<hull[0][2]<<"\n";
	cout<<oppositeSimplex[0]<<" "<<oppositeSimplex[1]<<" "<<oppositeSimplex[2];
	int k;
	cin>>k;
	auto weight2 = getweight(oppositeSimplex);
	auto projectionData2 = transformingConvexPolytopeForConvexDecomposition(inputData,oppositeSimplex);

	bool first = true;
	for(int i=0;i<projectionData1.first.size();i++){
		for(int j =0;j<projectionData1.first[i].size();j++){
			if(first){
				cout<<projectionData1.first[i][j];
				first=false;
			}
			else
				cout<<","<<projectionData1.first[i][j];
		}
		cout<<"\n";
		first = true;
	}
	first = true;
	cout<<"\n\n";
	for(int i=0;i<projectionData2.first.size();i++){
		for(int j =0;j<projectionData2.first[i].size();j++){
			if(first){
				cout<<projectionData2.first[i][j];
				first=false;
			}
			else
				cout<<","<<projectionData2.first[i][j];
		}
		cout<<"\n";
		first = true;

	}
	/*
	for(int i=0;i<hull.size();i++){
		for(int j =0;j<hull[i].size();j++){
			cout<<hull[i][j]<<" ";
		}
		cout<<"\n";
	}
	
	for(auto x: distanceMatrix){
		for(auto y :x){
			cout<<y<<" ";
		}
		cout<<"\n";
	}
	*/
	
}
