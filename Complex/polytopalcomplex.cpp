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


pair<vector<unsigned>,unsigned> findOppositeSimplex(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<vector<unsigned>> hull){
	vector<double> point1 = projectionData.second[0];
	vector<double> point2 = projectionData.second[1];
	double minweight1 = 9999;
	double minweight2 = 9999;
	unsigned projectionSimplex1,projectionSimplex2;
	int i=0;
	for(auto x : hull){
		vector<double> centroid(projectionData.first[0].size(),0.0);
		for(auto y :x){
			for(int j=0;j<projectionData.first[0].size();j++){
				centroid[j] +=((projectionData.first[y][j])/projectionData.first[0].size());
			}
		}
		auto d1 = utils::vectors_distance(centroid,point1);
		if(d1<minweight1){
			minweight1 = d1;
			projectionSimplex1 = i; 
		}
		auto d2 = utils::vectors_distance(centroid,point2);
		if(d2<minweight2){
			minweight2 = d2;
			projectionSimplex2 = i; 
		}
		i++;
	}
	if(minweight1<minweight2)
		return make_pair(hull[projectionSimplex1],0);
		
	return make_pair(hull[projectionSimplex2],1);
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
	pp1 = newpoints[newpoints.size()-2];
	pp2 = newpoints[newpoints.size()-1];
	auto trasformedData = utils::projectHSphere(inputData,center);
	pp.push_back(pp1);
	pp.push_back(pp2);
	pp.push_back(center);
	return make_pair(trasformedData,pp);
}
/*
def projectonplane(p1,p2,points1):
	cofficient = p1 - p2
	constant  = 0
	point_on_plane = p1
	for a,b in zip(point_on_plane,cofficient):
		constant = constant+a*b
	constant = (-1)*constant 
	
	pts1 = []
	for pt in points1:
		diff = pt-p2
		constterm =0;
		diffterm = 0;
		for x in range(0,len(cofficient)):
			constterm+= p2[x]*cofficient[x]
			diffterm+= diff[x]*cofficient[x]
		if((diffterm!=0).all()):
			t = ((-1)*constant - constterm)/diffterm
			pnt = t*(pt-p2)+p2
			pts1.append(pnt)
		else:
			pts1.append(pt)
	pts1 = np.array(pts1)
	return pts1
*/

vector<vector<double>> projectOnSimplexPlane(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<unsigned> hull,vector<double> ppoint){
	double scaledown = 100000;
	vector<double> pp;
	double d=utils::vectors_distance(ppoint,projectionData.second[2]);
	for(int i =0;i<ppoint.size();i++)
		pp.push_back(projectionData.second[2][i]+((ppoint[i]-projectionData.second[2][i])/(d*scaledown)));
	vector<double> normal;
	for(int i=0;i<projectionData.second[0].size();i++)
		normal.push_back(projectionData.second[0][i]-projectionData.second[1][i]);
	
	double constant  = 0;
	auto point_on_plane = projectionData.second[2];
	for (int i=0;i<point_on_plane.size();i++){
		constant += point_on_plane[i]*normal[i];
	}
	constant = (-1)*constant;
	
	vector<vector<double>> pts1;
	for (auto pt :projectionData.first){
		vector<double> diff;
		for(int i=0;i<pt.size();i++)
			diff.push_back(pt[i]-pp[i]);
		double constterm = 0;
		double diffterm = 0;
		for(int x=0;x<normal.size();x++){
			constterm+= pp[x]*normal[x];
			diffterm+= diff[x]*normal[x];
		}
		if(diffterm!=0){
			double t = ((-1)*constant - constterm)/diffterm;
			vector<double> pnt;
			for(int i=0;i<pt.size();i++){
				pnt.push_back(t*(pt[i]-pp[i])+pp[i]);
			}
			pts1.push_back(pnt);
		}
		else
			pts1.push_back(pt);
	}
	auto pts = utils::computePCA(pts1,pts1[0].size()-1);
	return pts.first;
}

polytopalComplex :: polytopalComplex(vector<vector<double>> &inputData){
	int k;
	populateDistanceMatrix(inputData);
	auto simplices = qdelaunay_o(inputData);
	auto hull = hullfromtriangulation(simplices);
	auto weight1 = getweight(hull[0]);
	auto projectionData1 = transformingConvexPolytopeForConvexDecomposition(inputData,hull[0]);
	auto oppositeSimplex = findOppositeSimplex(projectionData1,hull);
	vector<double> pp1,pp2;
	if(oppositeSimplex.second==0)
		pp1 = projectionData1.second[1];
	else
		pp1 = projectionData1.second[0];

	auto projectedonPlane1 = projectOnSimplexPlane(projectionData1,hull[0],pp1);
	//cout<<hull[0][0]<<" "<<hull[0][1]<<" "<<hull[0][2]<<"\n";
	//cout<<oppositeSimplex[0]<<" "<<oppositeSimplex[1]<<" "<<oppositeSimplex[2];
	
	//cin>>k;
	auto weight2 = getweight(hull[oppositeSimplex.second]);
	auto projectionData2 = transformingConvexPolytopeForConvexDecomposition(inputData,oppositeSimplex.first);
	if(utils::vectors_distance(pp1,projectionData2.second[0])<utils::vectors_distance(pp1,projectionData2.second[1]))
		pp2 = projectionData2.second[1];
	else
		pp2 = projectionData2.second[0];
	auto projectedonPlane2 = projectOnSimplexPlane(projectionData2,oppositeSimplex.first,pp2);

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
	first = true;
	cout<<"\n\n";
	for(int i=0;i<projectedonPlane1.size();i++){
		for(int j =0;j<projectedonPlane1[i].size();j++){
			if(first){
				cout<<projectedonPlane1[i][j];
				first=false;
			}
			else
				cout<<","<<projectedonPlane1[i][j];
		}
		cout<<"\n";
		first = true;

	}
	first = true;
	cout<<"\n\n";
	for(int i=0;i<projectedonPlane2.size();i++){
		for(int j =0;j<projectedonPlane2[i].size();j++){
			if(first){
				cout<<projectedonPlane2[i][j];
				first=false;
			}
			else
				cout<<","<<projectedonPlane2[i][j];
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
