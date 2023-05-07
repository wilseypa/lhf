#include<iostream>
#include "polytopalcomplex.hpp"

using namespace std;

double polytopalComplex :: distanceindex(unsigned x,unsigned y){
	if(x>y)
		return distanceMatrix[y][x-y];
	else
		return distanceMatrix[x][y-x];
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

std::vector<std::vector<unsigned>> polytopalComplex ::hullfromtriangulation(std::vector<std::vector<unsigned>> simplices){
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

pair<vector<vector<unsigned>>,vector<vector<unsigned>>> polytopalComplex ::neighbours(vector<unsigned>  polytop,std::vector<std::vector<unsigned>> simplices,int dim1,vector<vector<unsigned>> delaunaypart){   //report all the neighbouring polytopes that shares a face with polytope
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

pair<pair<vector<unsigned>,vector<vector<unsigned>>>,double> polytopalComplex ::mergeneighbors(vector<unsigned> polytop,std::vector<std::vector<unsigned>> simplices,vector<vector<double>>pointscoord,vector<vector<unsigned>> delaunaypart,set<unsigned> points){ //   merge neigbours simplices that results in maximum convex polytop with neighbors
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

pair<vector<unsigned>,set<unsigned>> polytopalComplex ::getnextsimplex(std::vector<std::vector<unsigned>> simplices,pair<vector<unsigned>,set<unsigned>> cSimpAddressedpoints,std::vector<std::vector<unsigned>> addressedsimplices){
	
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
pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> polytopalComplex ::optimalconvexization(std::vector<std::vector<unsigned>> simplices,vector<vector<double>> pointscoord,set<unsigned> originalpoints){   //Repeate the Maximum convexization for each simplex and returned the sorted list of convex polytopes by weight and vertices
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

pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> polytopalComplex :: sortofSizeandWeight(pair<pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>,pair<vector<int>,vector<double>>> df){

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
pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> polytopalComplex :: iterativeconvexization(std::vector<std::vector<unsigned>> simplices,int dim,vector<vector<double>> pointscoord){ //   Keep convex decomposition that minimizes convex polytopes required and minimizes maximum weight edge
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

	return make_pair(convexpartssorted,delaunayparts);
}
		

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> sortFaces(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexfaces){
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

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> pruneMaximalParts(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexFaces,pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexPolytopes){
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

pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  polytopalComplex :: informedConvexization(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  convexFaces, std::vector<std::vector<unsigned>> hull,pair<vector<vector<double>>,vector<double>> &projectionData){
	int d = projectionData.first[0].size()-1;
    pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> UPCF = convexFaces;
	int count=0;
	for(auto face : convexFaces.second){
		if(face.size()>2){
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
			auto sterographicProjection = projectOnSimplexPlane(ModifiedprojectionData,pp2,1);

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
pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>>  polytopalComplex ::generateConvexFaces(pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexFaces,std::vector<std::vector<unsigned>> hull, pair<vector<vector<double>>,vector<vector<double>>> &projectionData,vector<unsigned> projectionfacet,vector<double> pp){
	int d = projectionData.first[0].size()-1;
	auto sterographicProjection = projectOnSimplexPlane(projectionData,pp,1000);
	auto convexPolytopes = iterativeconvexization(hull,d,sterographicProjection);
	convexFaces = pruneMaximalParts(convexFaces,convexPolytopes);
	return convexFaces;
	
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

pair<vector<unsigned>,unsigned> polytopalComplex ::findOppositeSimplex(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<vector<unsigned>> hull){
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
pair<vector<vector<double>>,vector<double>> projectonCenteroidSphere(vector<vector<double>>inputData){
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


vector<vector<double>> polytopalComplex ::projectOnSimplexPlane(pair<vector<vector<double>>,vector<vector<double>>> projectionData,vector<double> ppoint,double scaledown){
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


	auto weight2 = getweight(hull[oppositeSimplex.second]);

	auto projectionData2 = transformingConvexPolytopeForConvexDecomposition(inputData,oppositeSimplex.first);

	if(utils::vectors_distance(pp1,projectionData2.second[0])<utils::vectors_distance(pp1,projectionData2.second[1]))
		pp2 = projectionData2.second[1];
	else
		pp2 = projectionData2.second[0];
		
	pair<vector<vector<unsigned>>,vector<vector<vector<unsigned>>>> convexFaces; 
	convexFaces = generateConvexFaces(convexFaces,hull,projectionData1,hull[0],pp1);
	convexFaces = generateConvexFaces(convexFaces,hull,projectionData2,oppositeSimplex.first,pp2);
	auto projectionData = projectonCenteroidSphere(inputData);
	convexFaces = informedConvexization(convexFaces, hull,projectionData);

	cout<<"Polytope Vertex Count\n\n";
	for(auto x: convexFaces.first)
		cout<<x.size()<<",";
	cout<<"\n\n";
	

}
