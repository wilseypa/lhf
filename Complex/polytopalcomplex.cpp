#include <iostream>
#include <numeric>
#include "polytopalcomplex.hpp"


using namespace std;

double polytopalComplex :: distanceindex(unsigned x,unsigned y){
	if(x>y)
		return distanceMatrix[y][x-y-1];
	else
		return distanceMatrix[x][y-x-1];
}
		
double polytopalComplex :: getweight(vector<unsigned> simplex){
	double weight =0;
	for(int i=0;i<simplex.size();i++){
		for(int j=i+1;j<simplex.size();j++){
			if(weight<distanceindex(simplex[i],simplex[j]))
				weight = distanceindex(simplex[i],simplex[j]);
		}
	}
	return weight;
}
void printData(vector<vector<double>> data){
    std::cout<<"\n\nData\n\n"<<flush;
	for(auto x: data){
		for(auto y:x){
			std::cout<<y<<","<<flush;
		}
	 std::cout<<"\n"<<flush;
	}
}
bool polytopalComplex :: checkCC_Simplex_Inclusion (std::vector<unsigned> simplex,std::vector<std::vector<double> >  &inputData,	std::vector<double> circumCenter){
	 int i = 0;
	 std::vector<std::vector<double>> matT;
	 std::vector<std::vector<double>> matPv;
	 for(int x =0 ;x<simplex.size()-1;x++){
			std::vector<double> tempmat;
			for(int j = 0;j<inputData[0].size();j++){
				tempmat.push_back(inputData[simplex[x]][j]-inputData[simplex[simplex.size()-1]][j]);
			}
			matT.push_back(tempmat);
	}
		for(int j = 0;j<inputData[0].size();j++){
				std::vector<double> tempmat;
				tempmat.push_back(circumCenter[j]-inputData[simplex[simplex.size()-1]][j]);
				matPv.push_back(tempmat);
			}
	 
	 std::vector<std::vector<double>> transposematT(matT.size(), std::vector<double> (matT[0].size(), 0));
	 for (int i = 0; i < matT.size(); ++i)
		for (int j = 0; j < matT[0].size(); ++j){
			transposematT[j][i]= matT[i][j];
		}
     std::vector<std::vector<double>> inverse = utils::inverseOfMatrix(transposematT,transposematT[0].size());
	 std::vector<std::vector<double>> lambda = utils::matrixMultiplication(inverse,matPv);
	 bool outside = false;
	 double sum = 0;
	 for(auto x:lambda){
		 sum+=x[0];
		 if(x[0]<0){
			outside = true;
			return outside;
		}
	 }
	 double lambda_n = 1-sum;
	 if(lambda_n<0)
		return true;
	 sum = sum + lambda_n;
	 if(sum >1)
		outside = true;
     return outside;
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
std::vector<std::vector<unsigned>> polytopalComplex ::generatesimplexfacets(std::vector<unsigned> c, unsigned k){
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

std::vector<std::vector<unsigned>> polytopalComplex ::hullfromtriangulation(std::vector<std::vector<unsigned>> &simplices){
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

vector<unsigned> polytopalComplex ::intersection(vector<unsigned> polytop,vector<unsigned> simplex){
 vector<unsigned> inter;
 for(auto y:simplex){
	 for(auto x:polytop){
		 if(y==x){
			inter.push_back(x);
			break;
		 }
	 }
 }
 return inter;
}

pair<vector<vector<unsigned>>,vector<vector<unsigned>>> polytopalComplex ::neighbours(vector<unsigned>  polytop,std::vector<std::vector<unsigned>> &simplices,int dim1,vector<vector<unsigned>> delaunaypart){   //report all the neighbouring polytopes that shares a face with polytope
	vector<vector<unsigned>> adjacency;
	for (auto y : simplices){
		sort(y.begin(),y.end());
		if(polytop != y){
			if((intersection(polytop,y)).size()==dim1)
				adjacency.push_back(y);
			if((intersection(polytop,y)).size()==dim1+1){
				bool present = false;
				for(auto x:delaunaypart)
					if(y ==x){
						present = true;
						break;
					}
					if(!present)
						delaunaypart.push_back(y);
				}
			}
	}
	return make_pair(adjacency,delaunaypart);
}

pair<pair<vector<unsigned>,vector<vector<unsigned>>>,double> polytopalComplex ::mergeneighborsPrev(vector<unsigned> polytop,std::vector<std::vector<unsigned>>& simplices,vector<vector<double>> &pointscoord,vector<vector<unsigned>> delaunaypart,set<unsigned> &points){ //   merge neigbours simplices that results in maximum convex polytop with neighbors
	auto neighborAndDelaunayParts = neighbours(polytop,simplices,pointscoord[0].size(),delaunaypart);
	auto neighbors = neighborAndDelaunayParts.first;
	delaunaypart = neighborAndDelaunayParts.second;	
	vector<unsigned> pts(points.begin(),points.end());
	double maxweightedge = 0;
	for(auto x : neighbors){
		set<unsigned> y(polytop.begin(),polytop.end());	
		int origSize = y.size();
		unsigned point;
		for(auto a:x){
			y.insert(a);
			if(y.size()>origSize){
				point = a;
				break;
			}
		}
		bool convex = true;
		auto hull = hullfromtriangulation(delaunaypart);
		for(auto d:hull){
			vector<vector<double>> coords;
			vector<unsigned> simplex;
			int k =0;
			for(auto i : d){
				auto it = find(pts.begin(), pts.end(), i);
				int index = it - pts.begin();
				coords.push_back(pointscoord[index]);
				simplex.push_back(k++);
			}
			simplex.push_back(k);
			auto it = find(pts.begin(), pts.end(), point);
			int index = it - pts.begin();
			coords.push_back(pointscoord[index]);
			for(auto p: polytop){
				if(find(d.begin(), d.end(), p) == d.end()){
					auto it = find(pts.begin(), pts.end(), p);
					int index = it - pts.begin();
					vector<double> pt = pointscoord[index];
					if(!checkCC_Simplex_Inclusion(simplex,coords,pt)){
						convex = false;
						break;	
					}
					if(!convex)
					break;	
				}
				if(!convex)
					break;	
			}
		 if(!convex)
		    break;		
		}
		if(convex){
			polytop.assign(y.begin(), y.end());
			if(maxweightedge<getweight(x))
				maxweightedge = getweight(x);
			delaunaypart.push_back(x);

		}
	}
	return make_pair(make_pair(polytop,delaunaypart),maxweightedge);
}
pair<pair<vector<unsigned>,vector<vector<unsigned>>>,double> polytopalComplex ::mergeneighbors(vector<unsigned> polytop,std::vector<std::vector<unsigned>>& simplices,vector<vector<double>> &pointscoord,vector<vector<unsigned>> delaunaypart,set<unsigned> &points){ //   merge neigbours simplices that results in maximum convex polytop with neighbors
	auto neighborAndDelaunayParts = neighbours(polytop,simplices,pointscoord[0].size(),delaunaypart);
	auto neighbors = neighborAndDelaunayParts.first;
	delaunaypart = neighborAndDelaunayParts.second;	
	vector<unsigned> pts(points.begin(),points.end());
	
	double maxweightedge = 0;
	for(auto x : neighbors){
		set<unsigned> y(polytop.begin(),polytop.end());
		for(auto a:x)
			y.insert(a);
		vector<vector<double>> coords;
		for(auto i : y){
			auto it = find(pts.begin(), pts.end(), i);
			int index = it - pts.begin();
			coords.push_back(pointscoord[index]);
		}
		auto simplices = qdelaunay_o(coords);
		auto hull = hullfromtriangulation(simplices);
		set<unsigned> points12;
		for(auto z : hull) 
			for(auto x : z){
				points12.insert(x);
			}
		if(points12.size()==y.size()){
			polytop.assign(y.begin(), y.end());
			if(maxweightedge<getweight(x))
				maxweightedge = getweight(x);
			delaunaypart.push_back(x);

		}
	}

	return make_pair(make_pair(polytop,delaunaypart),maxweightedge);
}

pair<vector<unsigned>,set<unsigned>> polytopalComplex ::getnextsimplex(std::vector<std::vector<unsigned>>& simplices,pair<vector<unsigned>,set<unsigned>> cSimpAddressedpoints,std::vector<std::vector<unsigned>> addressedsimplices){
	
	set<unsigned> addresspoints = cSimpAddressedpoints.second;
	for(auto a : cSimpAddressedpoints.first){
		addresspoints.insert(a);
	}
	vector<unsigned> addresspoint(addresspoints.begin(),addresspoints.end());
	std::vector<unsigned> x;
	bool present;
	for(auto y : simplices){
		x=y;
		present = false;
		for(auto k:addressedsimplices){
			sort(y.begin(),y.end());
			sort(k.begin(),k.end());
			if(y==k){
				present = true;
				break;
			}
		}
		if(!present)
			break;
	}
	return make_pair(x,addresspoints);
}
pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> polytopalComplex ::optimalconvexization(std::vector<std::vector<unsigned>> &simplices,vector<vector<double>> &pointscoord,set<unsigned> &originalpoints){   //Repeate the Maximum convexization for each simplex and returned the sorted list of convex polytopes by weight and vertices
	vector<vector<unsigned>> convexparts;
	vector<double> maxdist;
	vector<vector<vector<unsigned>>> delaunayparts;
	vector<int> sizecd;
	set<unsigned> totalpoints;
	int k;
	vector<unsigned> origpts(originalpoints.begin(),originalpoints.end());
	auto updatedsimplices = simplices;
	pair<vector<unsigned>,set<unsigned>> cSimpAddressedpoints;
	vector<unsigned> x;
	std::vector<std::vector<unsigned>> addressedsimplices;
	
	for(int i =0;i<simplices.size();i++){
		vector<unsigned> psimplex = cSimpAddressedpoints.first;
		cSimpAddressedpoints.first = x;
		cSimpAddressedpoints = getnextsimplex(simplices,cSimpAddressedpoints,addressedsimplices);
		
		if(psimplex==cSimpAddressedpoints.first)
			break;
		x = cSimpAddressedpoints.first;
		sort(x.begin(), x.end());
		vector<vector<unsigned>> delaunaypart;
		delaunaypart.push_back(x);
		double distance =0;
		double maxweight = getweight(x);
		double maxval = 0;
		while(true){
			auto x1delaunaypartw = mergeneighbors(x,simplices,pointscoord,delaunaypart,originalpoints);
			if(maxweight<x1delaunaypartw.second)
				maxweight = x1delaunaypartw.second;
			set<unsigned> xx;
			set<unsigned> pp;
			for(auto p:x1delaunaypartw.first.second)
				addressedsimplices.push_back(p);
			for(auto p:x1delaunaypartw.first.first){
				xx.insert(p);
			}
			for(auto p:x){
				pp.insert(p);
			}
			if(xx==pp)
				break;
			x = x1delaunaypartw.first.first;
			delaunaypart = x1delaunaypartw.first.second;
			
		}
		bool present = false;
		for(auto y :convexparts)
			if(x==y){
				present = true;
				break;
			}
	
		if(!present){
			double maxval=0;
			double distance = 0;
			convexparts.push_back(x);
			delaunayparts.push_back(delaunaypart);
			/*
			vector<vector<double>> coords;
			
			for(auto i : x){
				auto it = find(origpts.begin(), origpts.end(), i);
				int index = it - origpts.begin();
				coords.push_back(pointscoord[index]);

			}
			auto simplices = qdelaunay_o(coords);
			auto hull = hullfromtriangulation(simplices);

			for(auto p : hull){
				auto it = find(origpts.begin(), origpts.end(), x[p[0]]);
				int index1 = it - origpts.begin();
				it = find(origpts.begin(), origpts.end(), x[p[1]]);
				int index2 = it - origpts.begin();
				double d = utils::vectors_distance(pointscoord[index1],pointscoord[index2]);
				distance = distance + d;
				if(maxval < d)
					maxval = d;
			}
			*/
			maxdist.push_back(maxval);
			sizecd.push_back(x.size());
		}
		set<vector<unsigned>> k;
		for (auto tp : delaunayparts){
			for(auto y:tp){
				sort(y.begin(),y.end());
				k.insert(y);
			}
		}
		if(k.size()==simplices.size())
			break;
	}
	
	
	return make_pair(make_pair(convexparts,delaunayparts),make_pair(sizecd,maxdist));
}
bool sortDecending(const pair<int,int> &a,const pair<int,int> &b)
{
    return (a.first > b.first);
}

pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> polytopalComplex :: sortofSizeandWeight(pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> &df){

	vector<pair<double,unsigned>>	sortbyWeight;
	unsigned i =0;
	for(double k : df.second.second){
		sortbyWeight.push_back(make_pair(k,i));
		i++;
	}
	sort(sortbyWeight.begin(), sortbyWeight.end());
	vector<pair<unsigned,unsigned>>	sortbysize;
	for(auto k : sortbyWeight){
		sortbysize.push_back(make_pair(df.second.first[k.second],k.second));
	}
	sort(sortbysize.begin(), sortbysize.end(), sortDecending);
	
	vector<vector<unsigned>> convexparts;
	vector<double> maxdist;
	vector<vector<vector<unsigned>>> delaunayparts;
	vector<int> sizecd;
	
	for(auto k : sortbysize){
		convexparts.push_back(df.first.first[k.second]);
		delaunayparts.push_back(df.first.second[k.second]);
		sizecd.push_back(df.second.first[k.second]);
		maxdist.push_back(df.second.second[k.second]);
	}
	
	
	return make_pair(make_pair(convexparts,delaunayparts),make_pair(sizecd,maxdist));
}
pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex :: iterativeconvexization(std::vector<std::vector<unsigned>>& simplices,int dim,vector<vector<double>> &pointscoord){ //   Keep convex decomposition that minimizes convex polytopes required and minimizes maximum weight edge
	vector<int> valid;
	vector<vector<unsigned>> convexpartsunion;
	vector<vector<vector<unsigned>>> delaunaypartsfinal;
	set<unsigned> originalpoints;
	for (auto k :simplices){
		for(auto i : k){
			originalpoints.insert(i);
		}
	}
	vector<double> weights;
	vector<vector<vector<unsigned>>> delaunayparts;
	auto remainingsimplices = simplices;
	int premaining = 0;
	while(true){
		auto df = optimalconvexization(remainingsimplices,pointscoord,originalpoints);
		
		df = sortofSizeandWeight(df);
	
		unsigned i=0;
		vector<unsigned> pointsaddressed;
		for(unsigned i=0;i<df.first.first.size();i++){
			auto x = df.first.first[i];
			auto y = df.first.second[i];
			auto weight = df.second.second[i];
			bool check = true;
			for(auto t:df.first.second[i]){
					for( auto tt : delaunaypartsfinal){
						for(auto p:tt){
							sort(t.begin(),t.end());
							sort(p.begin(),p.end());
							if(t == p)
								check = false;
							}
					if(!check)	
						break;
				}
				if(!check)
					break;
			}
			if(check){
				valid.push_back(i);
				convexpartsunion.push_back(x);
				delaunayparts.push_back(y);
				weights.push_back(weight);
				delaunaypartsfinal.push_back(y);
			}
		}
		remainingsimplices.clear();
		int remaining = 0;
		for(auto i : simplices){
			bool present = false;
			for(auto x : delaunaypartsfinal){
				for(auto y:x){
					sort(y.begin(),y.end());
					sort(i.begin(),i.end());
					if (y==i){
						present = true;
						break;
					}
				if(present)
					break;	
				}
				if(present)
					break;	
			}
			if(!present){
				remaining = remaining+1;
				remainingsimplices.push_back(i);
			}
		}	
		if(remaining==premaining){
			if(remaining !=0)
				for (auto r : remainingsimplices){
					convexpartsunion.push_back(r);
					weights.push_back(getweight(r));
					vector<vector<unsigned>> rr;
					rr.push_back(r);
					delaunayparts.push_back(rr);
					delaunaypartsfinal.push_back(rr);
				}
			break;
		}
		premaining = remaining;
	}
	vector<vector<unsigned>>  convexpartssorted;
	for(auto x : convexpartsunion){
		sort(x.begin(),x.end());
		convexpartssorted.push_back(x);
	}
/*
	cout<<"\n";
	for(auto data : pointscoord){
		for(auto coord:data)
			cout<<coord<<",";
		cout<<"0\n";
	}
	for(auto y : delaunayparts){
		for(auto x : y){
			bool t = true;
			for(auto yy:x)
				if(t){
					std::cout<<yy<<flush;
					t = false;
				}
				else
					std::cout<<","<<yy<<flush;
			std::cout<<"\n";
			}
		std::cout<<"0,0,0\n"<<flush;
	}
	int kk;
	cin>>kk;
*/
	return make_pair(convexpartssorted,delaunayparts);
}
		

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex ::sortFaces(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexfaces){
	vector<vector<unsigned>> newconvexfaces;
	vector<vector<vector<unsigned>>> newconvexfacesdel;
	
	vector<pair<unsigned,unsigned>>	sortbysize;
	unsigned i =0;
	for(auto k : convexfaces.first){
		sortbysize.push_back(make_pair(k.size(),i));
		i++;
	}
	sort(sortbysize.begin(), sortbysize.end(), sortDecending);
	for(auto x : sortbysize){
		newconvexfaces.push_back(convexfaces.first[x.second]);
		newconvexfacesdel.push_back(convexfaces.second[x.second]);
	}
	return make_pair(newconvexfaces,newconvexfacesdel);
}

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex ::pruneMaximalParts(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexFaces,pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexPolytopes){
	vector<vector<unsigned>> newconvexF1 = convexFaces.first;
	vector<vector<vector<unsigned>>> newconvexfd1 = convexFaces.second;
	for(unsigned i=0;i< convexPolytopes.first.size();i++){
		newconvexF1.push_back(convexPolytopes.first[i]);
		newconvexfd1.push_back(convexPolytopes.second[i]);
	}
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexFaces1 = make_pair(newconvexF1,newconvexfd1);
	convexFaces1 =  sortFaces(convexFaces1);
	vector<vector<unsigned>> newconvexF;
	vector<vector<vector<unsigned>>> newconvexfd;
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> newconvexFaces = make_pair(newconvexF,newconvexfd);

	for(unsigned i=0;i< convexFaces1.first.size();i++){
		bool present = false;
		for(auto x : convexFaces1.second[i]){
			for(unsigned j =0;j<newconvexFaces.first.size();j++){
				for(auto y : newconvexFaces.second[j]){
					sort(x.begin(),x.end());
					sort(y.begin(),y.end());
					if(x==y){
						present = true;
						break;
					}
				}
				if(present)
					break;
			}
			if(present)
				break;
		}
		if(!present){
			newconvexFaces.first.push_back(convexFaces1.first[i]);
			newconvexFaces.second.push_back(convexFaces1.second[i]);
		}
	}

	return newconvexFaces;
}

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  polytopalComplex :: informedConvexization(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  &convexFaces, std::vector<std::vector<unsigned>> &hull,pair<vector<vector<double>>,vector<double>> &projectionData){
	int d = projectionData.first[0].size()-1;
    pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> UPCF = convexFaces;
	int count=0;
	for(auto face : convexFaces.second){
		if(face.size()>=1){
			vector<double> centeroid(d+1,0.0);
			for( auto b : convexFaces.first[count]){
				for(int i=0;i<d+1;i++){
					centeroid[i] +=(projectionData.first[b][i]/convexFaces.first[i].size());
				}
			}

			vector<vector<double>> pp;	
			pp.push_back(centeroid);
			pp = utils::projectHSphere(pp,projectionData.second);	

			vector<double> pp2(projectionData.first[0].size(),0.0);
			for(int i=0;i<projectionData.first[0].size();i++){
				pp2[i] = 2*projectionData.second[i] - pp[0][i];
			}

			pp.push_back(pp2);
			pp.push_back(projectionData.second);		
			pair<vector<vector<double>>,vector<vector<double>>> ModifiedprojectionData;
			ModifiedprojectionData.first = projectionData.first;
			ModifiedprojectionData.second =pp;
			auto sterographicProjection = projectOnSimplexPlane(ModifiedprojectionData,pp2,100000);
			auto convexPolytopes = iterativeconvexization(hull,d,sterographicProjection);
			UPCF = pruneMaximalParts(UPCF,convexPolytopes);

		}
		count++;
	}
	for(auto x : hull){
		bool present = false;
		for (auto y : UPCF.second){
			for(auto z : y){
				sort(x.begin(),x.end());
				sort(z.begin(),z.end());
				if(x==z){
					present = true;
					break;
				}
			}
			if(present)
				break;
		}
		if(!present){
			vector<vector<unsigned>> parts;
			parts.push_back(x);
			UPCF.first.push_back(x);
			UPCF.second.push_back(parts);

		}
	}
	
	return UPCF;
}
pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  polytopalComplex ::generateConvexFaces(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> &convexFaces,std::vector<std::vector<unsigned>> hull, pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<unsigned> projectionfacet,vector<double> pp){
	int d = projectionData.first[0].size()-1;
	auto sterographicProjection = projectOnSimplexPlane(projectionData,pp,100000);
	auto convexPolytopes = iterativeconvexization(hull,d,sterographicProjection);
	return 	pruneMaximalParts(convexFaces,convexPolytopes);
}
vector<vector<double>> polytopalComplex ::projectpointsOnUnitdSphere(vector<vector<double>> &inputData,vector<double> &center){
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

pair<vector<unsigned>,unsigned> polytopalComplex ::findOppositeSimplex(pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<vector<unsigned>> hull){
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
pair<vector<vector<double>>,vector<double>> polytopalComplex :: projectonCenteroidSphere(vector<vector<double>>&inputData){
	vector<double> centeroid(inputData[0].size(),0.0);
	for( auto b : inputData){
		for(int i=0;i<inputData[0].size();i++){
			centeroid[i] +=(b[i]/inputData.size());
		}
	}
	auto trasformedData = utils::projectHSphere(inputData,centeroid);
	return make_pair(trasformedData,centeroid);
}

pair<vector<vector<double>>,vector<vector<double>>> polytopalComplex ::transformingConvexPolytopeForConvexDecomposition(vector<vector<double>> &inputData,vector<unsigned> projectionSimplex){
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


vector<vector<double>> polytopalComplex ::projectOnSimplexPlane(pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<double> ppoint,double scaledown){
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
std::vector<std::vector<unsigned>> polytopalComplex ::  transformhull(std::vector<std::vector<unsigned>> originalhull,set<unsigned> polytopeIndices){ // Global to Local Indices
	std::vector<std::vector<unsigned>> newhull = originalhull;
	int index=0;
	for(auto x: polytopeIndices){
		for(int i=0;i<newhull.size();i++){
			for(int j=0;j<newhull[i].size();j++)
				if(newhull[i][j] == x){
					newhull[i][j] = index;
				}
		}
		index++;
	}
	return newhull;
}
pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex ::  transformCF(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexFaces,set<unsigned> polytopeIndices){ // local to global Indices
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> newcf= convexFaces;
	int index=polytopeIndices.size()-1;
	for(auto x= polytopeIndices.rbegin();x!=polytopeIndices.rend();x++){
		for(int i=0;i<newcf.first.size();i++){
			for(int j=0;j<newcf.first[i].size();j++)
				if(newcf.first[i][j] == index){
					newcf.first[i][j] = *x;
				}
			
			for(int k=0;k<newcf.second[i].size();k++)
				for(int l=0;l<newcf.second[i][k].size();l++)
					if(newcf.second[i][k][l] == index){
						newcf.second[i][k][l] = *x;
					}
		}
		index--;
	}
	return newcf;
}

double polytopalComplex :: getmaxWeight(vector<vector<unsigned>> &delparts){
	double maxWeight = 0;
	for(auto x:delparts)
		if(maxWeight < getweight(x))	
			maxWeight = getweight(x);
	return maxWeight;
}

double average(std::vector<double> const& v){
    if(v.empty()){
        return 0;
    }

    auto const count = static_cast<double>(v.size());
    return std::reduce(v.begin(), v.end()) / count;
}
template <typename T>
dd_MatrixPtr polytopalComplex :: dd_PolyFile2Matrix_2(std::vector<std::vector<T>> A, std::vector<T> B, dd_ErrorType *Error)
{
  (*Error) = dd_NoError;
  dd_MatrixPtr M = NULL;
  long m_input, i;
  long d_input, j;
  dd_RepresentationType rep;
  mytype value;
  int newformat = dd_FALSE, successful = dd_FALSE, linearity = dd_FALSE;
  char comsave[dd_linelenmax], numbtype[dd_wordlenmax];
  dd_NumberType NT;

  dd_init(value);
  rep = dd_Inequality;
  newformat = dd_TRUE;
  linearity = dd_TRUE;

  // numbtype from typedata
  strcpy(numbtype, "integer");
  NT = dd_GetNumberType(numbtype);
  if (NT == dd_Unknown)
  {
    (*Error) = dd_ImproperInputFormat;
    goto _L99;
  }
  m_input = A.size();
  d_input = A[0].size() + 1;
  M = dd_CreateMatrix(m_input, d_input);
  M->representation = rep;
  M->numbtype = NT;

  for (i = 0; i < m_input; i++)
  {
    dd_set_d(value, B[i]);
    dd_set(M->matrix[i][0], value);
    for (j = 1; j < d_input; j++)
    {
      dd_set_d(value, A[i][j - 1]);
      dd_set(M->matrix[i][j], value);
    }
  }
  //dd_WriteMatrix(stdout, M);
  successful = dd_TRUE;
  strcpy(comsave, "linearity 1 3");
  
_L99:;
  dd_clear(value);
  return M;
}

vector<vector<double>> polytopalComplex :: HyperplaneToVertexRepresentation(std::vector<std::vector<double>> &AA,std::vector<double>& BB){
  std::vector<std::vector<double>> v;
  dd_ErrorType err = dd_NoError;
  dd_MatrixPtr M = NULL;
  dd_set_global_constants();

  /* Read data from stdin */
  M = dd_PolyFile2Matrix_2(AA, BB, &err);

  dd_PolyhedraPtr poly;
  dd_MatrixPtr A;
  /* compute the second representation */
  poly = dd_DDMatrix2Poly(M, &err);

  switch (poly->representation)
  {
  case dd_Inequality:
    A = dd_CopyGenerators(poly);
    for (int i = 0; i < A->rowsize; i++)
    {
      std::vector<double> temp;
      for (int j = 0; j < A->colsize; j++)
        temp.push_back(*(A->matrix[i][j]));
      v.push_back(temp);
    }
    // dd_WriteMatrix(stdout,A);
    dd_FreeMatrix(A);
    break;

  case dd_Generator:
    //fprintf(stdout, "The second representation:\n");
    A = dd_CopyInequalities(poly);
    dd_WriteMatrix(stdout, A);
    dd_FreeMatrix(A);
    break;

  default:
    break;
  }
    std::vector<std::vector<double>> vv;
	for (int i = 0; i < v.size(); i++)
	{   std::vector<double> v1;
		for (int j = 1; j < v[0].size(); j++)
			v1.push_back(-1*v[i][j]);
		vv.push_back(v1);
	}
  return vv;
}
pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex :: postprocess(std::vector<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>> toprocess){
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> processedFaces;
	vector<vector<unsigned>> polytopes;
	vector<vector<vector<unsigned>>> polydelparts;
	for(auto x: toprocess){
		set<unsigned> indices;
		vector<vector<unsigned>> delparts;
		for(auto xx : x.first){
			for(auto ind :xx)
				indices.insert(ind);
		}
		for(auto xx : x.second){
			for(auto delpa :xx)
				delparts.push_back(delpa);
		}
		vector<unsigned> inde(indices.begin(),indices.end());
		polytopes.push_back(inde);
		polydelparts.push_back(delparts);
	}
	
	return make_pair(polytopes,polydelparts);
}


pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex ::  flattened(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> faces,vector<vector<double>> coords,int dim){
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> flatten;
    vector<vector<unsigned>> minimal_facets;
    
    
	for(auto x: faces.first){
		for(auto y:x){
				cout<<y<<",";
		}
		cout<<"\n";
	}
	
    set<unsigned> vertices;	
    int i=0,j=0;
	for( auto x :faces.first){
		set<unsigned> boundary;	
		cout<<"\nnextface::\n";		
		for(auto xx : x){
			vector<int> count;
			cout<<"\nxx::"<<xx;
			for(int j=i+1; j <faces.first.size();j++){
				for(auto yk : faces.second[j])
					for(auto yy : yk){
						if(xx==yy){
							count.push_back(j);
					}
				}
			}
			cout<<count.size()<<"ROhit\n";
			if(count.size() ==1){
				boundary.insert(xx);
				vertices.insert(xx);
			}
			if(count.size()>=2){
			   set<int> unique;
			   for(auto x:count)
					unique.insert(x);
				if(unique.size()>=2)
				{
					boundary.insert(xx);
					vertices.insert(xx);
				}
			}
		}
		vector<unsigned> bound(boundary.begin(),boundary.end());
		minimal_facets.push_back(bound);
	}
	vector<unsigned> vert(vertices.begin(),vertices.end());
	vector<vector<double>> data;
	for(auto x: vert)
		data.push_back(coords[x]);
    auto simplices = qdelaunay_o(data);
	set<unsigned> hull;
	for(auto x:simplices){
		for(auto k:x){
			hull.insert(vert[k]);
		}
	}
	vector<unsigned> hullbound(hull.begin(),hull.end());
	minimal_facets.push_back(hullbound);

	for(auto x: minimal_facets){
		for(auto y:x){
			cout<<y<<" ";
		}
		cout<<"\n";
	}
	int k;
	//cin>>k;
	return flatten;
}


pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex :: approxDecomposition(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> faces,vector<vector<double>> coords,int dim, double approxscale){
	int k;
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> unprocessedFaces=sortFaces(faces);
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> approxFaces;

	if(faces.first.size()<=1)
		return faces;
	int i=0;
	while(true){
	i=0;
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> processedFaces;
	vector<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>> approxFaces1;
	while(unprocessedFaces.first.size()>0){
		auto x = unprocessedFaces.first[0];
		pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> Neighbors;
		int j=0;
		for(auto y:unprocessedFaces.first){
			if(x!=y){
			   if(intersection(x,y).size()>=dim){
			      Neighbors.first.push_back(y);
			      Neighbors.second.push_back(unprocessedFaces.second[j]);
			   }
			}
			j++;
		}
		auto mergeResults = mergeApprox(make_pair(x,unprocessedFaces.second[0]),Neighbors,approxscale,coords);
		processedFaces.first.push_back(x);
		processedFaces.second.push_back(unprocessedFaces.second[0]);
		int k=0;
		for(auto p: mergeResults.first){
			processedFaces.first.push_back(p);
			processedFaces.second.push_back(mergeResults.second[k]);
			k++;
		}
		mergeResults.first.push_back(x);
		mergeResults.second.push_back(unprocessedFaces.second[0]); 
		approxFaces1.push_back(mergeResults);
	    i++;
	    unprocessedFaces = updateUnprocessed(unprocessedFaces,processedFaces);
	    approxFaces = postprocess(approxFaces1);

	}
	if(processedFaces.first.size() == approxFaces.first.size()){
	   std::cout<<"Reduction:: "<<approxFaces.first.size()<<"\n";
	   break;
   }
	else{
		unprocessedFaces = sortFaces(approxFaces);
	}
	}
	
	std::cout<<"ConvexFaces::"<<"\n";
	i=0;
	for(auto x:approxFaces.first){
		int iii=0;
		cout<<"[";
		for(auto y: approxFaces.second[i]){
			cout<<"[";
			int ii =0;
			for(auto k:y){
				if(ii==0)
					std::cout<<k;
				else
					std::cout<<","<<k;
				ii =1;
			}
			if(iii==approxFaces.second[i].size()-1)
				std::cout<<"]";
			else
				std::cout<<"],";
			iii =iii+1;
		}
		std::cout<<"]\n";
		i++;
	}
	std::cout<<"\n";
	//std::cin>>k;
	
	
	//return flattened(approxFaces,coords,dim);
	return approxFaces;

}

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex ::updateUnprocessed(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> unprocessed, pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> processed){
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> updatedunprocessed;
	
	for(int i=0;i<unprocessed.first.size();i++){
		bool present = false;
		for(auto k=0;k<processed.first.size();k++){
		    if(unprocessed.first[i]==processed.first[k]){
				present = true;
				break;
			}
		}
		if(!present){
			updatedunprocessed.first.push_back(unprocessed.first[i]);
			updatedunprocessed.second.push_back(unprocessed.second[i]);
		}
	}
	return updatedunprocessed;
}

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex :: mergeApprox(pair<vector<unsigned>,vector<vector<unsigned>>> x,pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> Neighbors,double approxscale,vector<vector<double>> coords){
	
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> approxMerged;
	for(auto y:x.first)
		cout<<y<<",";
	
	cout<<"\n Neighbors \n";
	for(auto k:Neighbors.first){
		for(auto t:k){
			cout<<t<<",";
		}
		cout<<"\n";
	}
	cout<<"\nMerge \n";
	int i=0;
	for(auto k:Neighbors.first){
		auto yy = intersection(x.first,k);
		bool merge = true;
		std::cout<<"\n"<<yy.size()<<" "<<k.size()<<" Intersection"<<"\n";
		for(auto t:k){
			    bool atleastone = false;
				for(auto y:yy){
				std::cout<<"\nDistance : "<<utils::vectors_distance(coords[y],coords[t])<<"\n Polytop\n";
				for(auto l:coords[y]){
					cout<<l<<",";
				}
				cout<<"\n Neighbor";
				for(auto l:coords[t]){
					cout<<l<<",";
				}
				cout<<"\n";
				
				if(utils::vectors_distance(coords[y],coords[t])<approxscale){
				    atleastone = true;
				    break;
				}
			}
			if(!atleastone){
				merge = false;
				break;
			}
		}
		if(merge){
		std::cout<<"\nValid::";
		  approxMerged.first.push_back(k);
		  approxMerged.second.push_back(Neighbors.second[i]);
		  for(auto kk:k)
			cout<<kk<<",";
		}
		cout<<"\n";
		i++;
	}
	
	int k;
	//cin>>k;
	return approxMerged;
}


polytopalComplex :: polytopalComplex(vector<vector<double>> &inputData,double approxfactor){
   /*std::vector<std::vector<double>> A = {{1.235, 0.235, 0.3246}, {0.235, 1.125, 0.678}, {0.124124, 0.1252346, 1.657}, {-1.01, -1.231, -1.234523}, {-1.3426, -1.125, 0.658}, {0.568, -1.12321, -1.235423}, {-1.2365, 0.345423, -1.568}, {-1.124, -1.34562, -1.24}};
   std::vector<double> B = {0.26, 0.1243, 0.346, 1.124, 1.123, 1.6758, 1.123, 2.325};
   auto v = HyperplaneToVertexRepresentation(A,B);

	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < v[0].size(); j++)
			std::cout<<v[i][j]<<" ";
		std::cout<<std::endl;
	}
	int k;
	std::cout<<std::endl;
    std::cin>>k;
    */
	auto simplices = qdelaunay_o(inputData);
	std::cout<<"Delaunay  Size ::"<<simplices.size()<<"\n";
	auto mdata = inverseStereoGraphicProjection(inputData);
	populateDistanceMatrix(inputData);
	dim = inputData[0].size()+1;
	dataSize =inputData.size();
	/*
	std::cout<<"Simplice Size ::"<<simplices.size()<<"\n";
	for(auto y : simplices){
		bool t = true;
			for(auto yy:y)
				if(t){
					std::cout<<yy<<flush;
					t = false;
				}
				else
					std::cout<<","<<yy<<flush;
		std::cout<<"\n0,0,0\n"<<flush;
	}
	*/
	auto convexPolytopes1 = iterativeconvexization(simplices,inputData[0].size(),inputData);
	std::cout<<"convexPoly  Size ::"<<convexPolytopes1.second.size()<<"\n";
	int i=0;
	/*
	for(auto poly :convexPolytopes1.second){
		std::vector<std::vector<double>> facetData;
		for(auto x : convexPolytopes1.first[i]){
			std::cout<<x<<" ";
		}
		std::cout<<"\n\n\n";
		for(auto x : convexPolytopes1.first[i]){
			facetData.push_back(mdata[x]);
			for(auto kk : mdata[x])
				std::cout<<kk<<",";
			std::cout<<"\n";
		}
		std::cout<<"\n\n\n";
		for(auto y : poly){
			bool t = true;
			for(auto yy:y)
				if(t){
					std::cout<<yy<<flush;
					t = false;
				}
				else
					std::cout<<","<<yy<<flush;
			std::cout<<"\n"<<flush;
		}
		auto HP = utils::computePCA(facetData,3);
		std::cout<<"\n\n\n";
		for(auto x : HP.first){
			for(auto kk : x)
				std::cout<<kk<<",";
			std::cout<<"\n";
		}
		std::cout<<"\n\n\n";
		for(auto x : HP.second){
			for(auto kk : x)
				std::cout<<kk<<",";
			std::cout<<"\n";
		}
		int kkk;
		std::cin>>kkk;
		std::cout<<"\n0,0,0\n"<<flush;
		i++;
	}
	*/
	
	set<polytope,cmp> polys;
	polytope temp;
	set<unsigned> p;
	for(int i=0;i<inputData.size();i++)
		p.insert(i);
	temp.polytopeIndices = p;
	temp.coordinates = mdata;
	temp.Deltriangulation = simplices;
	polys.insert(temp);
	i =0;
	int level = 0;
	polytopalArrayList.push_back(polys);
	while(true){
		set<polytope,cmp> polys;
		unsigned cofaceI=0;
		cout<<"\nConvex Parts Size  :: Level ::"<<level<<" has "<<polytopalArrayList[level].size()<<" Size \n"<<flush;
		for(auto poly:polytopalArrayList[level]){
		//for(auto poly:polytopalArrayList[level]){
			std::cout<<poly.polytopeIndices.size()<<","<<flush;
			if(poly.polytopeIndices.size()>poly.coordinates[0].size()+1){
				auto originalhull = poly.Deltriangulation; //Original Indices
				auto hull = transformhull(originalhull,poly.polytopeIndices); // local Indices
				auto weight1 = getweight(hull[0]);
				auto projectionData1 = transformingConvexPolytopeForConvexDecomposition(poly.coordinates,hull[0]);
				auto oppositeSimplex = findOppositeSimplex(projectionData1,hull);
				vector<double> pp1,pp2;
				if(oppositeSimplex.second==0)
					pp1 = projectionData1.second[1];
				else
					pp1 = projectionData1.second[0];
				auto weight2 = getweight(hull[oppositeSimplex.second]);
				auto projectionData2 = transformingConvexPolytopeForConvexDecomposition(poly.coordinates,oppositeSimplex.first);
				if(utils::vectors_distance(pp1,projectionData2.second[0])<utils::vectors_distance(pp1,projectionData2.second[1]))
					pp2 = projectionData2.second[1];
				else
					pp2 = projectionData2.second[0];
				pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexFaces;
				convexFaces = generateConvexFaces(convexFaces,hull,projectionData1,hull[0],pp1);
				convexFaces = generateConvexFaces(convexFaces,hull,projectionData2,oppositeSimplex.first,pp2);
		
				auto projectionData = projectonCenteroidSphere(poly.coordinates);
				convexFaces = informedConvexization(convexFaces, hull,projectionData);
				
				if(level ==0)
					convexFaces= pruneMaximalParts(convexFaces,convexPolytopes1);
				
				convexFaces = transformCF(convexFaces,poly.polytopeIndices); // local Indices
				int i=0;
				convexFaces = approxDecomposition(convexFaces,inputData,dim-level-1,approxfactor);
				for(auto poly :convexFaces.first){
					polytope temp;
					set<unsigned> p(poly.begin(),poly.end());
					temp.polytopeIndices = p;
					vector<vector<double>> coords;
					for(auto x:temp.polytopeIndices)
						coords.push_back(mdata[x]);
					auto pca = utils:: computePCA(coords,coords[0].size());
					temp.coordinates = pca.first;
					std::cout<<"\nEigen Vectors::\n";
					for(auto xy : pca.second){
						for(auto kk : xy)
							std::cout<<kk<<",";
						std::cout<<"\n";
					}
					std::cout<<"\nMean::\n";
					for(auto xxx: utils::transpose(coords))
						std::cout<<average(xxx)<<",";
					std::cout<<"\n";
					for(auto x :temp.polytopeIndices)
					   std::cout<<x<<" ";
					  std::cout<<"\n";
					for(auto y : temp.coordinates){
						bool t = true;
						for(auto yy : y)
						if(t){
							 t = false;
							std::cout<<yy<<flush;
							}
						else
							std::cout<<yy<<","<<flush;
						cout<<"\n";
					}
					cout<<"0,0,0\n";
					temp.Deltriangulation = hullfromtriangulation(convexFaces.second[i]);
					cout<<"\n";
					for(auto y : hullfromtriangulation(convexFaces.second[i])){
						bool t = true;
						for(auto yy : y)
						 if(t){
							 t = false;
							std::cout<<yy<<flush;
							}
						else
							std::cout<<","<<yy<<flush;
						cout<<"\n";
					}
					cout<<"0,0,0\n";

					temp.weight = getmaxWeight(convexFaces.second[i]);
					temp.cofaceIndices.insert(cofaceI);
					auto pos = polys.find(temp);
					if (pos == polys.end())
						polys.insert(temp);
					else{
						auto it = next(polys.begin(), distance(polys.begin(), pos));
						auto polyupdate = *it;
						polyupdate.cofaceIndices.insert(cofaceI);
						polys.erase(it);
						polys.insert(polyupdate);
					}
					i++;
				}
			}
			else{
				//Do for simplices
				int i=0;
				vector<unsigned> pp(poly.polytopeIndices.begin(),poly.polytopeIndices.end());
				for(auto poly :generatesimplexfacets(pp,poly.polytopeIndices.size()-1)){
					polytope temp;
					set<unsigned> p(poly.begin(),poly.end());
					temp.polytopeIndices = p;
					temp.weight = getweight(poly);
					temp.cofaceIndices.insert(cofaceI);
					auto pos = polys.find(temp);
					if (pos == polys.end())
						polys.insert(temp);
					else{
						auto it = next(polys.begin(), distance(polys.begin(), pos));
						auto polyupdate = *it;
						polyupdate.cofaceIndices.insert(cofaceI);
						polys.erase(it);
						polys.insert(polyupdate);
					}
					i++;
				}
			}
			cofaceI++;
		}
		polytopalArrayList.push_back(polys);
		level++;
		if(level+1 == dim-1)
			break;
	}
	set<polytope,cmp> polysE;
	unsigned cofaceI =0;
	for(auto poly:polytopalArrayList[level]){
		for(auto p : poly.Deltriangulation){
			polytope temp;
			set<unsigned> pp(p.begin(),p.end());
			temp.polytopeIndices = pp;
			temp.Deltriangulation.push_back(p);
			temp.weight = getweight(p);
			temp.cofaceIndices.insert(cofaceI);
			auto pos = polysE.find(temp);
			if (pos == polysE.end())
				polysE.insert(temp);
			else{
				 auto it = next(polysE.begin(), distance(polysE.begin(), pos));
				 auto polyupdate = *it;
				 polyupdate.cofaceIndices.insert(cofaceI);
				 polysE.erase(it);
				 polysE.insert(polyupdate);
			}
		}
		cofaceI++;
	}
	polytopalArrayList.push_back(polysE);
	level++;
    set<polytope,cmp> polysV;
	cofaceI =0;
	for(auto poly:polytopalArrayList[level]){
		for(auto p : poly.polytopeIndices){
			polytope temp;
			temp.polytopeIndices.insert(p);
			temp.cofaceIndices.insert(cofaceI);
			temp.weight = 0.0;
					temp.cofaceIndices.insert(cofaceI);
			auto pos = polysV.find(temp);
			if (pos == polysV.end())
				polysV.insert(temp);
			else{
				 auto it = next(polysV.begin(), distance(polysV.begin(), pos));
				 auto polyupdate = *it;
				 polyupdate.cofaceIndices.insert(cofaceI);
				 polysV.erase(it);
				 polysV.insert(polyupdate);
			}
		}
		cofaceI++;
	}
	polytopalArrayList.push_back(polysV);
	int l=0;
	for(auto xLevel:polytopalArrayList){
		std::cout<<"\nLevel "<<l++<<"has Size "<<xLevel.size()<<"\n["<<flush;
		for(auto poly :xLevel){
			std::cout<<"[["<<flush;
			for(auto y : poly.polytopeIndices)
				std::cout<<y<<","<<flush;
			std::cout<<"]"<<flush;
			std::cout<<"["<<flush;
			for(auto y : poly.cofaceIndices)
				std::cout<<y<<","<<flush;
			std::cout<<"]["<<flush;
			for(auto y : poly.Deltriangulation){
				cout<<"(";
				for(auto yy : y)
					std::cout<<yy<<","<<flush;
				}
			std::cout<<")"<<poly.weight<<"]]->"<<flush;
		}	
	}
	int row=0;
	int col=0;
	/*
	for(auto x: distanceMatrix){
		std::cout<<"\n";
		col = 0;
		for(auto y:x){
			std::cout<<"("<<row<<","<<col<<")"<<y<<" ";
			col++;
			}
		row++;
		}
	*/
}
std::vector<polytope> polytopalComplex :: getCofacetModified(set<unsigned> cofaceIndices,int codim){
	std::set<polytope> cofaces;
	auto cofacesdim = polytopalArrayList[codim];
	int cofaceI =0;
	for(auto x:cofaceIndices){
		auto it = next(cofacesdim.begin(), x);
		auto polytop = (*it);
		for(auto del : polytop.Deltriangulation){
			polytope temp;
			set<unsigned> pp(del.begin(),del.end());
			temp.polytopeIndices = pp;
			temp.Deltriangulation.push_back(del);
			temp.weight = getweight(del);
			temp.cofaceIndices.insert(cofaceI);
			auto pos = cofaces.find(temp);
			if (pos == cofaces.end())
				cofaces.insert(temp);
			else{
				 auto it = next(cofaces.begin(), distance(cofaces.begin(), pos));
				 auto polyupdate = *it;
				 polyupdate.cofaceIndices.insert(cofaceI);
				 cofaces.erase(it);
				 cofaces.insert(polyupdate);
			}		
		}
		cofaceI++;
	}
	vector<polytope> cof(cofaces.begin(),cofaces.end());
	return cof	;
}

std::vector<polytope> polytopalComplex :: getCofacet(set<unsigned> cofaceIndices,int codim){
	std::vector<polytope> cofaces;
	auto cofacesdim = polytopalArrayList[codim];
	for(auto x:cofaceIndices){
		auto it = next(cofacesdim.begin(), x);
		cofaces.push_back(*it);
	}
	return cofaces;
}

set<polytope>  polytopalComplex :: generateFaces(int d){
	set<polytope> faces;
	int cofaceI = 0;
	for(auto poly = polytopalArrayList[d].rbegin(); poly != polytopalArrayList[d].rend(); poly++){
		auto polytop = (*poly);
		for(auto del : polytop.Deltriangulation){
			polytope temp;
			set<unsigned> pp(del.begin(),del.end());
			temp.polytopeIndices = pp;
			temp.Deltriangulation.push_back(del);
			temp.weight = getweight(del);
			temp.cofaceIndices.insert(cofaceI);
			auto pos = faces.find(temp);
			if (pos == faces.end())
				faces.insert(temp);
			else{
				 auto it = next(faces.begin(), distance(faces.begin(), pos));
				 auto polyupdate = *it;
				 polyupdate.cofaceIndices.insert(cofaceI);
				 faces.erase(it);
				 faces.insert(polyupdate);
			}		
		}
		cofaceI++;
	}	
	return faces;

}


void polytopalComplex :: persistenceMatrix(){
	
	for(auto level = 4;level>=0;level--){
	     cout<<"\n"<<polytopalArrayList[level].size()<<"\n";
	}
	int level = 4;
	for(int i=0;i<3;i++){
		vector<vector<int>> matrix(polytopalArrayList[level-i].size(), vector<int>(polytopalArrayList[level-1-i].size(),0));
		int k=polytopalArrayList[level-i].size()-1;
		for(auto x : polytopalArrayList[level-i]){
			for(auto y :x.cofaceIndices)
				matrix[k][y] = 1;
			k--;	
		}
		for(int k =0; k<polytopalArrayList[level-1-i].size();k++){
repeat:
			int pivot = 0;
			for(int c = 0;c<matrix.size();c++)
				if(matrix[c][k]==1)
				   pivot = c;
			for(int kk = k-1; kk>=0;kk--){
				bool reduce = false;
				if(matrix[pivot][kk]==1){
					reduce = true;
					for(int p = pivot+1;p<matrix.size();p++)
					    if(matrix[p][kk]==1){
							reduce =false;
							break;
						}
				}
				if(reduce){
					for(int pp = 0;pp<matrix.size();pp++){
						matrix[pp][k] = (matrix[pp][k] + matrix[pp][kk])%2;
					}
				goto repeat;
				}
			}
		}
		cout<<"here"<<flush;
		int pr=0,pc=0;
		for(pc = 0;pc<matrix[0].size();pc++){
			int row = -1;
			for(pr =0; pr<matrix.size();pr++)
				if(matrix[pr][pc]==1)
					row = pr;
			if(row != -1){
				 auto it1 = next(polytopalArrayList[level-i].begin(), matrix.size()-row-1);
				 auto it2 = next(polytopalArrayList[level-i-1].begin(), matrix[0].size()-pc-1);
				if((*it2).weight != (*it1).weight){
					cout<<"\nbetties "<<(*it1).weight <<" "<<(*it2).weight;
				}
			}	 
		}
		
	}
}


/*
std::vector<polytope> polytopalComplex :: persistenceByDim(std::vector<polytope> pivots,int d){
	typename std::vector<polytope>::iterator it = pivots.begin();

	std::vector<polytope> nextPivots;	 	//Pivots for the next dimension
	std::unordered_map<polytope, std::vector<polytope>,MyHashFunction> v;	//Store only the reduction matrix V and compute R implicity
	std::unordered_map<polytope, polytope,MyHashFunction> pivotPairs;	//For each pivot, which column has that pivot
	//Iterate over columns to reduce in reverse order
	for(auto columnIndexIter = polytopalArrayList[dim-d].rbegin(); columnIndexIter != polytopalArrayList[dim-d].rend(); columnIndexIter++){

		polytope poly = (*columnIndexIter);	//The current simplex

		//Not a pivot -> need to reduce
		if(it == pivots.end() || (*it).weight != poly.weight || (*it).polytopeIndices != poly.polytopeIndices){
			std::vector<polytope> faceList = getCofacetModified(poly.cofaceIndices,dim-d-1);
			std::vector<polytope> columnV;	
			columnV.push_back(poly); 
			std::make_heap(faceList.begin(), faceList.end());

			while(true){
				polytope pivot;
				while(!faceList.empty()){		
					pivot = faceList.front();
					//Rotate the heap
					std::pop_heap(faceList.begin(), faceList.end());
					faceList.pop_back();
					if(!faceList.empty() && pivot.polytopeIndices == faceList.front().polytopeIndices){ //Coface is in twice -> evaluates to 0 mod 2
						//Rotate the heap
						std::pop_heap(faceList.begin(), faceList.end());
						faceList.pop_back();
					} else{
						faceList.push_back(pivot);
						std::push_heap(faceList.begin(), faceList.end());
						break;
					}		
				}
				if(faceList.empty()){ //Column completely reduced
					break;
				} else if(pivotPairs.find(pivot) == pivotPairs.end()){ //Column cannot be reduced
					pivotPairs.insert({pivot, poly});
					nextPivots.push_back(pivot);		
					std::sort(columnV.begin(), columnV.end());
					auto it = columnV.begin();
					while(it != columnV.end()){
						if((it+1) != columnV.end() && (*it)==*(it+1)) ++it;
						else v[poly].push_back(*it);
						++it;
					}

					if(poly.weight != pivot.weight){
						vector<double> des = {d, std::min(pivot.weight, poly.weight), std::max(pivot.weight, poly.weight)};
						bettiTable.push_back(des);
					}

					break;
				} else{ 
			
					for(polytope p : v[pivotPairs[pivot]]){
						columnV.push_back(p);
						std::vector<polytope> faces = getCofacetModified(p.cofaceIndices,dim-d-1);
						faceList.insert(faceList.end(), faces.begin(), faces.end());
					}
					std::make_heap(faceList.begin(), faceList.end());
				}
			}
		} else ++it;
	}

	return nextPivots;
}
*/
vector<vector<double>> polytopalComplex :: rescaledataandcenteraroundorigin(vector<vector<double>> &data){
	double maxvalue = 99999;
	vector<vector<double>> newdata;
	vector<double> maxim(data[0].size(),0.0);
	vector<double> minim(data[0].size(),maxvalue);
	for(auto y :data){
		for(int j=0;j<data[0].size();j++){
			if(maxim[j] < y[j]){
				maxim[j] = y[j];
			}
			if(minim[j] > y[j]){
				minim[j] = y[j];
			}
		}
	}
	for(auto y :data){
		vector<double> temp;
		for(int j=0;j<data[0].size();j++){
			temp.push_back(((y[j] -minim[j])/(maxim[j]-minim[j])) -0.5);
		}
		newdata.push_back(temp);
	}
	
	return newdata;
}
vector<vector<double>> polytopalComplex :: inverseStereoGraphicProjection(vector<vector<double>> &data){
	vector<vector<double>> newdata;
	auto updateddata = rescaledataandcenteraroundorigin(data);
	double radius = 0.1;
	vector<double> origin(data[0].size(),0.0);
	vector<double> center(data[0].size(),0.0);
	vector<double> projection(data[0].size(),0.0);
    center.push_back(radius);
    projection.push_back(2*radius);
    origin.push_back(0);
		
	for(auto pt : data){
		pt.push_back(0);
		double scaler = (0.1/utils::vectors_distance(center,pt));
		vector<double> temp;
		for(int j=0;j<=data[0].size();j++){
			double coord = (center[j])*(1-scaler) + pt[j]*scaler;
			temp.push_back(coord);
			cout<<coord<<" ";
		}
		cout<<"\n";
		newdata.push_back(temp);
	}
    return newdata;
}
/*
void polytopalComplex :: persistence(){
	//Get all edges for the simplexArrayList or simplexTree
	auto edges = polytopalArrayList[dim-1];

	if(edges.size() <= 1) return;
	std::unordered_map<unsigned, unsigned> mappedIndices;	
	std::vector<polytope> pivots; 
	unsigned mstSize = 0;
	unionFind uf(dataSize);
	for(auto edgeIter = edges.rbegin(); edgeIter != edges.rend(); edgeIter++){
		std::set<unsigned>::iterator it = (*edgeIter).polytopeIndices.begin();
		if( mappedIndices.size() == 0 || mappedIndices.find(*it) == mappedIndices.end() ) mappedIndices.insert( std::make_pair(*it, mappedIndices.size()) );
		int c1 = uf.find(mappedIndices.find(*it)->second);
		it++;
		if( mappedIndices.find(*it) == mappedIndices.end() ) mappedIndices.insert( std::make_pair(*it, mappedIndices.size()) );
		int c2 = uf.find(mappedIndices.find(*it)->second);
		if(c1 != c2){
			uf.join(c1, c2);
			mstSize++;
			pivots.push_back((*edgeIter));
			vector<double> des = { 0, 0, (*edgeIter).weight};
			bettiTable.push_back(des);
		}
		if(mstSize >= edges.size()-1) break;
	}

	for(int i=0; i<dataSize; i++){
		if(uf.find(i) == i){ //i is the name of a connected component
			vector<double> des = { 0, 0, 1000};
			bettiTable.push_back(des);
		}
	}
	for(unsigned d = 1; d < dim && d < polytopalArrayList.size()-1; d++){
		pivots = persistenceByDim(pivots, d);
	}
	return;
}

*/

std::vector<polytope> polytopalComplex :: persistenceByDim(std::vector<polytope> pivots,int d){
	typename std::vector<polytope>::iterator it = pivots.begin();

	std::vector<polytope> nextPivots;	 	//Pivots for the next dimension
	std::unordered_map<polytope, std::vector<polytope>,MyHashFunction> v;	//Store only the reduction matrix V and compute R implicity
	std::unordered_map<polytope, polytope,MyHashFunction> pivotPairs;	//For each pivot, which column has that pivot
	//Iterate over columns to reduce in reverse order
	for(auto columnIndexIter = polytopalArrayList[dim-d].rbegin(); columnIndexIter != polytopalArrayList[dim-d].rend(); columnIndexIter++){

		polytope poly = (*columnIndexIter);	//The current simplex

		//Not a pivot -> need to reduce
		if(it == pivots.end() || (*it).weight != poly.weight || (*it).polytopeIndices != poly.polytopeIndices){
			std::vector<polytope> faceList = getCofacet(poly.cofaceIndices,dim-d-1);
			std::vector<polytope> columnV;	
			columnV.push_back(poly); 
			std::make_heap(faceList.begin(), faceList.end());

			while(true){
				polytope pivot;
				while(!faceList.empty()){		
					pivot = faceList.front();
					//Rotate the heap
					std::pop_heap(faceList.begin(), faceList.end());
					faceList.pop_back();
					if(!faceList.empty() && pivot.polytopeIndices == faceList.front().polytopeIndices){ //Coface is in twice -> evaluates to 0 mod 2
						//Rotate the heap
						std::pop_heap(faceList.begin(), faceList.end());
						faceList.pop_back();
					} else{
						faceList.push_back(pivot);
						std::push_heap(faceList.begin(), faceList.end());
						break;
					}		
				}
				if(faceList.empty()){ //Column completely reduced
					break;
				} else if(pivotPairs.find(pivot) == pivotPairs.end()){ //Column cannot be reduced
					pivotPairs.insert({pivot, poly});
					nextPivots.push_back(pivot);		
					std::sort(columnV.begin(), columnV.end());
					auto it = columnV.begin();
					while(it != columnV.end()){
						if((it+1) != columnV.end() && (*it)==*(it+1)) ++it;
						else v[poly].push_back(*it);
						++it;
					}

					if(poly.weight != pivot.weight){
						vector<double> des = {d, std::min(pivot.weight, poly.weight), std::max(pivot.weight, poly.weight)};
						bettiTable.push_back(des);
					}

					break;
				} else{ 
			
					for(polytope p : v[pivotPairs[pivot]]){
						columnV.push_back(p);
						std::vector<polytope> faces = getCofacet(p.cofaceIndices,dim-d-1);
						faceList.insert(faceList.end(), faces.begin(), faces.end());
					}
					std::make_heap(faceList.begin(), faceList.end());
				}
			}
		} else ++it;
	}

	return nextPivots;
}

void polytopalComplex :: persistence(){
	//Get all edges for the simplexArrayList or simplexTree
	auto edges = polytopalArrayList[dim-1];

	if(edges.size() <= 1) return;
	std::unordered_map<unsigned, unsigned> mappedIndices;	
	std::vector<polytope> pivots; 
	unsigned mstSize = 0;
	unionFind uf(dataSize);
	for(auto edgeIter = edges.rbegin(); edgeIter != edges.rend(); edgeIter++){
		std::set<unsigned>::iterator it = (*edgeIter).polytopeIndices.begin();
		if( mappedIndices.size() == 0 || mappedIndices.find(*it) == mappedIndices.end() ) mappedIndices.insert( std::make_pair(*it, mappedIndices.size()) );
		int c1 = uf.find(mappedIndices.find(*it)->second);
		it++;
		if( mappedIndices.find(*it) == mappedIndices.end() ) mappedIndices.insert( std::make_pair(*it, mappedIndices.size()) );
		int c2 = uf.find(mappedIndices.find(*it)->second);
		if(c1 != c2){
			uf.join(c1, c2);
			mstSize++;
			pivots.push_back((*edgeIter));
			vector<double> des = { 0, 0, (*edgeIter).weight};
			bettiTable.push_back(des);
		}
		if(mstSize >= edges.size()-1) break;
	}

	for(int i=0; i<dataSize; i++){
		if(uf.find(i) == i){ //i is the name of a connected component
			vector<double> des = { 0, 0, 1000};
			bettiTable.push_back(des);
		}
	}
	for(unsigned d = 1; d < dim && d < polytopalArrayList.size()-1; d++){
		pivots = persistenceByDim(pivots, d);
	}
	return;
}


