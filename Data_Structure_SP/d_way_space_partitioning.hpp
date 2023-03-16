
/* 
 * Need for data structure for d-way space partioning:
 * Issues with current partioning trees:
 * 1. Existing data structures have worst time complexity in higher dimensions.
 * 2. Data points required explods with dimension 2^d.
 * 3. Have same time complexity as bruit force in higher dimensions.
 * 4. Partioning of space increases with dimensions 2^d
 * 
 * How d-way trees can be fruitful in higher dimensionsl space.
 * 1. The space partioning is linear to the data dimension.
 * 2. Required point is O(d+1) as compared to O(2^d).
 * 3. Is a way to provide standary multi-way partioning.
 * 4. NN search complexity can be decreased.
 */ 


/*
				*************************Polytopal Complex***************************
*1. Will be useful in convex partioning of space.
*2. Could be an alternative/assist in faster approximate Delaunay Triangulation.
*3. Partion space in convex components.
*4. Highly Scalable with NN searches.

*/
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iterator>
#include <bits/stdc++.h>
#include <string>
#include "../Utils/utils.hpp"
#include "../Utils/readInput.hpp"


#include <chrono>
using namespace std::chrono;
    
#define MAXRADIUS 9999

double precision = 10000;
std::vector<std::vector<double>> referenceHypertetrhedron;
int dim;
std::vector<double> time4;
std::vector<double> time1;
std::vector<double> time2;
std::vector<double> time3;
void makeCombiUtil(std::vector<std::vector<int> >& ans,std::vector<int>& tmp, int n, int left, int k);
std::vector<std::vector<int> > makeCombi(int n, int k);

struct lexical_compare_points {
	bool operator() (const std::vector<double> a, const std::vector<double> b) const {
		for(int i=0;i<a.size();i++)
		    if(a[i]!=b[i])
		         return a[i] > b[i];
		return false;    
	}
};

struct comp_by_radius{
	bool operator() (const std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>> a, const std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>> b) const {
		return  a.second.second < b.second.second;   
	}
};
class dwaytreenode {
public:
    std::vector<double> coordinates;
    std::vector<std::vector<double>> directionVectors;
    double radius;
    dwaytreenode* parent;
    int parenttochilddirection;
    std::vector<dwaytreenode*> children;
    dwaytreenode()
    {
		this->coordinates.assign(dim, 0);
		this->parent = nullptr;
    }
    dwaytreenode(std::vector<double> coords,double radius,int direction)
    {
        this->coordinates = coords;
        this->radius = radius;
        this->parenttochilddirection = direction;
    }
    dwaytreenode(std::vector<double> coords,double radius,int direction,std::vector<std::vector<double>> directionVectors)
    {
        this->coordinates = coords;
        this->radius = radius;
        this->parenttochilddirection = direction;
        this->directionVectors = directionVectors;
    }
    dwaytreenode* buildDwayTree(std::vector<std::vector<double>> data,int direction,std::string,int);
    void printTree(dwaytreenode *root);
    void printLevelOrder(dwaytreenode*,dwaytreenode* root);
    void printCurrentLevel(dwaytreenode*,dwaytreenode* root, int level);
    int height(dwaytreenode* node);
    void pathToCell(std::vector<dwaytreenode*> &path,dwaytreenode* root,dwaytreenode* node);
	std::vector<std::vector<double>>  findCellBoundingPolytope(dwaytreenode* root, dwaytreenode* node);
	dwaytreenode*  findNearestNeighbor(dwaytreenode* root, std::vector<double>);
	bool checkPointInBall(dwaytreenode* root, std::vector<double>,double,std::vector<std::vector<double>>);
	std::set<std::vector<double>> pointInBall(dwaytreenode* root, std::vector<double>,double);
	std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> meshGeneration(dwaytreenode* mainroot,dwaytreenode* root, double beta, int homologydim,double epsilon);
	std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> mergedmesh(dwaytreenode* root, std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> meshestomerge, double beta, int homologydim,double epsilon);
	std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> filterValidSimplices(dwaytreenode* root,std::set<std::set<std::vector<double>,lexical_compare_points>> simplicestocheck,double beta);
	std::set<std::set<std::vector<double>,lexical_compare_points>> generateNewSimplices(std::set<std::set<std::vector<double>,lexical_compare_points>> simplices, std::set<std::vector<double>,lexical_compare_points> points);
	std::set<std::vector<double>,lexical_compare_points> generatePoints(std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> partition);
	std::set<std::set<std::vector<double>,lexical_compare_points>> generateCombinations(std::set<std::vector<double>,lexical_compare_points> isolatedpoints,double homologydim);
	std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> generateAllSimplicestoCheck(std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> meshestoMerge,int homologydim);
	std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> validatesimplices(dwaytreenode* root,std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> tovalidate,double beta);
	std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> generatePartitionFromIsolatedPoints(dwaytreenode* root,std::set<std::vector<double>,lexical_compare_points> isolatedpoints,int homologydim,double epsilon,double beta);
	std::pair<bool,std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>> validatesimplex(dwaytreenode* root,std::set<std::vector<double>,lexical_compare_points> simplex,int homologydim,double);
	std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> generateAllSimplicestoCheck1(dwaytreenode* root,std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> propermeshestomerge,std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> newpartition,int homologydim,double);
	std::set<std::set<std::vector<double>,lexical_compare_points>> generateCombinationsall(std::set<std::set<std::vector<double>,lexical_compare_points>> allsimplices,double homologydim);
	std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> stichMesh(dwaytreenode* mainroot,dwaytreenode* root, std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> meshestomerge,double beta, int homologydim,double epsilon);
	std::vector<std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>>> computeSpaceFabricNeighbourhood(dwaytreenode* root,std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> ,double epsilon,int);
	std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> generateSimplicesToConsider(dwaytreenode* root,std::vector<std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>>> stichNeighboorhood,std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> propermeshestomerge,std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> newpartition,int homologydim,double epsilon,double beta);
	std::set<std::set<std::vector<double>,lexical_compare_points>> generateKSimplices(std::vector<std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>>> stichNeighboorhood,int k);
	std::set<std::set<std::vector<double>,lexical_compare_points>> generateKSimplices(std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> stichNeighboorhood,int k,int homologydim);
	bool checkInsertSubDsimplex(std::set<unsigned> dsimplex,std::vector<std::vector<double>> data,std::vector<std::vector<double>> distMatrix,double beta,dwaytreenode* tree,std::string tp);

};

void dwaytreenode:: printLevelOrder(dwaytreenode* originalroot,dwaytreenode* root){
    int h = height(root);
    int i;
    for (i = 1; i <= h; i++){
//		std::cout<<"Level ::"<<i<<"\n";
        printCurrentLevel(originalroot,root, i);
	}
}
 
void dwaytreenode:: printCurrentLevel(dwaytreenode* originalroot,dwaytreenode* root, int level){
    if (root == nullptr)
        return;
    
    if (level == 1){
		for(auto xx: root->coordinates){
			    std::cout<<xx<<" ";
		}
		std::cout<<root->radius<<"\n";
		return;
		/*
        auto coord = originalroot->findCellBoundingPolytope(originalroot,root);
        std::vector<std::vector<double>> coords;
        int enume = 0;
        */ 
     /*   for(auto xx: root->coordinates){
			    std::cout<<xx<<" ";
		}
		std::cout<<"\n";*/
        /*for(auto x:coord){
			if(x.size()>0)
			   coords.push_back(x);
			else{
			   
			    int xxx=0;
			   std::vector<double> temp;
			   for(auto xx: root->coordinates){
			    temp.push_back(referenceHypertetrhedron[enume][xxx]+xx);
			    xxx++;
				}
				
			   coords.push_back(temp);
			}
			enume++;   
		}
		*/
/*		for(auto tpy:coords){
			for(auto t : tpy){
				std::cout<<t<<" ";
			}
			std::cout<<"\n";
		}
	*/
		/*	
        int count = 0;
        int k;
    	for(auto x : coords){
			/*std::cout<<"Begin\n";
			for(auto t : x){
				std::cout<<t<<" ";
			}
			std::cout<<"End\n";
			*/
			/*
			if(x.size()>0){
			std::vector<std::vector<double>> pts;
			pts.push_back(x);
			int cnt = 0;
			for(auto axis:referenceHypertetrhedron){
				if(cnt!=count){
					int g = 0;
					std::vector<double> updatedaxis;
					for(auto shift:x){
						double num = double(shift)+double(axis[g]);
						updatedaxis.push_back(num);
						g++;
						}
					pts.push_back(updatedaxis);
				}
				cnt++;
			}
			auto axisall = makeCombi(referenceHypertetrhedron.size()-1,referenceHypertetrhedron.size()-2);
			for(auto comb : axisall){
				std::vector<std::vector<double>> hyplanvertices;
				hyplanvertices.push_back(x);
				for(auto y : comb){
					hyplanvertices.push_back(pts[y]);
				}
				* */
		/*	for(auto tt: hyplanvertices){
				for(auto t : tt){
					std::cout<<t<<" ";
				}
			std::cout<<"End\n";
			}
			*/
			/*
			auto hyperplane = utils::generateHyperplaneFromVertices(hyplanvertices,root->coordinates);
			* */
	/*		for(auto y : hyperplane.first)
				  std::cout<<y<<",";
			std::cout<<hyperplane.second<<"\n";
	*/
	//		}	
	//	}
	//	count++;
//	}
//	std::cin>>k;
//		std::cout<<"\n************************\n";
	}
    else if (level > 1) {
		for(auto child:root->children){
			printCurrentLevel(originalroot,child, level - 1);

		}
    }
}

int dwaytreenode:: height(dwaytreenode* node){
    if (node == nullptr)
        return 0;
    else {
		int maxheight = 0;
		for(auto child:node->children){
			int height1 = height(child);
			if(height1>maxheight)
				maxheight = height1;
        }
            return (maxheight + 1);
    }
}


std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> dwaytreenode::filterValidSimplices(dwaytreenode* root,std::set<std::set<std::vector<double>,lexical_compare_points>> simplicestocheck,double beta){
	std::set<std::set<std::vector<double>,lexical_compare_points>> validlist;
	
	std::cout<<simplicestocheck.size()<<"\n";
	std::set<unsigned> simplex;
	if(simplicestocheck.size()>0)
		for(int i=0;i<(*simplicestocheck.begin()).size();i++)
			simplex.insert(i);
	std::set<std::vector<double>,lexical_compare_points> pointsAccounted;
	std::set<std::vector<double>,lexical_compare_points> totalPoints;
	for(auto x:simplicestocheck){
		std::vector<std::vector<double>> simp(x.begin(), x.end());
		std::vector<std::vector<double>> distMatrix(x.size(),std::vector<double>(x.size()));
		int i=0,j=0;
		for(auto y:x){
			j=0;
			for(auto z:x){
				distMatrix[i][j] = utils::vectors_distance(y,z);
				j++;
			}
			i++;
		}
		
		auto cc =  utils::circumCenter(simplex,simp);
		auto radius = sqrt(utils::circumRadius(simplex,&distMatrix));
		if(!checkPointInBall(root, cc,radius, simp)){
			validlist.insert(x);
			for(auto pt:x){
				pointsAccounted.insert(pt);
			}
		}
		
		for(auto pt:x){
				totalPoints.insert(pt);
			}
			
	}
	std::vector<std::vector<double>> tp(totalPoints.begin(),totalPoints.end());

	for(auto x : pointsAccounted)
		tp.erase(std::remove(tp.begin(), tp.end(), x), tp.end());
	std::set<std::vector<double>,lexical_compare_points> isolatedpoints(tp.begin(),tp.end());
	
	//Need to filter out list based on beta inclusion rule.
	//Need to report isolated points also
	return std::make_pair(validlist,isolatedpoints);
}

std::set<std::set<std::vector<double>,lexical_compare_points>> dwaytreenode:: generateNewSimplices(std::set<std::set<std::vector<double>,lexical_compare_points>> simplices,std::set<std::vector<double>,lexical_compare_points> points){
	
	std::set<std::set<std::vector<double>,lexical_compare_points>> newsimplices;
	
	for(auto x: simplices){
		auto facets = generateCombinations(x,x.size()-2);
		for(auto z : facets){
			for(auto y : points){
				z.insert(y);
				newsimplices.insert(z);
				z.erase(y);
			}
		}
	}
	
	//It will generate simplices by taking facets from A with each points in B
	return newsimplices;
}

std::set<std::vector<double>,lexical_compare_points> dwaytreenode:: generatePoints(std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> partition){
	std::set<std::vector<double>,lexical_compare_points> points;
	
	for(auto x: partition.second)
		points.insert(x);
	
	for(auto y: partition.first)
		for(auto z:y)
			points.insert(z);
	
	//It will degenerate partition to its points
	return points;
}

std::set<std::set<std::vector<double>,lexical_compare_points>> dwaytreenode:: generateCombinations(std::set<std::vector<double>,lexical_compare_points> isolatedpoints,double homologydim){
    std::set<std::set<std::vector<double>,lexical_compare_points>> simplices;
    
    auto simp = makeCombi(isolatedpoints.size(), homologydim+1);
    for(auto x : simp){
		std::set<std::vector<double>,lexical_compare_points> sim;
		for(auto y:x){
			sim.insert(*std::next(isolatedpoints.begin(), y));
		}
		simplices.insert(sim);
	}
    //generate simplices from isolate points.
    return simplices;
}

std::set<std::set<std::vector<double>,lexical_compare_points>> dwaytreenode:: generateCombinationsall(std::set<std::set<std::vector<double>,lexical_compare_points>> allsimplices,double homologydim){
    std::set<std::set<std::vector<double>,lexical_compare_points>> simplices;
    for(auto isolatedpoints : allsimplices){
		auto simp = makeCombi(isolatedpoints.size(), homologydim+1);
		for(auto x : simp){
			std::set<std::vector<double>,lexical_compare_points> sim;
			for(auto y:x){
				sim.insert(*std::next(isolatedpoints.begin(), y));
			}
			simplices.insert(sim);
		}
	}
    //generate simplices from isolate points.
    return simplices;
}


std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> dwaytreenode:: generateAllSimplicestoCheck(std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> meshestoMerge, int homologydim){

std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> simplices;

for(int k = 0;k<homologydim;k++){
   int otherk = homologydim-k-1;
   std::set<std::set<std::vector<double>,lexical_compare_points>> facetsa;
   std::set<std::set<std::vector<double>,lexical_compare_points>> facetsb;
   if(k==0){
	   for(auto v:meshestoMerge[0].second){
			std::set<std::vector<double>,lexical_compare_points> faa;
			faa.insert(v);
			facetsa.insert(faa);
		}
   }
   if(otherk==0){
	   for(auto v:meshestoMerge[1].second){
		   	std::set<std::vector<double>,lexical_compare_points> fbb;
			fbb.insert(v);
			facetsb.insert(fbb);
		}
   }
   for(auto sima: meshestoMerge[0].first){
        auto facets = generateCombinations(sima,k);
        for(auto f:facets)
            facetsa.insert(f);
   }
   for(auto simb: meshestoMerge[1].first){
        auto facets = generateCombinations(simb,otherk);
		for(auto f:facets)
            facetsb.insert(f);
   }
   std::set<unsigned> repsimplex;
   for(int i=0;i<homologydim+1;i++)
			repsimplex.insert(i);
			
   std::set<std::set<std::vector<double>,lexical_compare_points>> finalsimplices;   
   for(auto firsthalf:facetsa){
	   for(auto secondhalf:facetsb){
		 std::set<std::vector<double>,lexical_compare_points> simplex;
		 std::set_union(firsthalf.begin(), firsthalf.end(),secondhalf.begin(), secondhalf.end(),std::inserter(simplex, simplex.begin()));
         std::vector<std::vector<double>> simp(simplex.begin(), simplex.end());
		 std::vector<std::vector<double>> distMatrix(simplex.size(),std::vector<double>(simplex.size()));
		 int i=0,j=0;
		 for(auto y:simplex){
			j=0;
			for(auto z:simplex){
				distMatrix[i][j] = utils::vectors_distance(y,z);
				j++;
			}
			i++;
		 } 
		
		 auto cc =  utils::circumCenter(repsimplex,simp);
		 auto radius = sqrt(utils::circumRadius(repsimplex,&distMatrix));
		 simplices.insert(std::make_pair(simplex,std::make_pair(cc,radius)));
	   }
   }
}

return simplices;

}

std::pair<bool,std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>> dwaytreenode:: validatesimplex(dwaytreenode* root,std::set<std::vector<double>,lexical_compare_points> simplex,int homologydim,double beta){
	
	std::vector<std::vector<double>> simpl(simplex.begin(), simplex.end());
	std::vector<std::vector<double>> distMatrix(simplex.size(),std::vector<double>(simplex.size()));
	std::set<unsigned> repsimplex;
	for(int i=0;i<homologydim+1;i++)
		repsimplex.insert(i);
	int i= 0;
	int j= 0;
	for(auto y:simpl){
		j=0;
		for(auto z:simpl){
			distMatrix[i][j] = utils::vectors_distance(y,z);
			j++;
		}
		i++;
	}
	auto start1 = high_resolution_clock::now();
	auto cc =  utils::circumCenter(repsimplex,simpl);
	auto stop1 = high_resolution_clock::now();
    auto duration1 = duration_cast<microseconds>(stop1 - start1);
    time1.push_back(duration1.count());
   
	//auto radius = sqrt(utils::circumRadius(repsimplex,&distMatrix));
	auto radius = utils::vectors_distance(simpl[0],cc);
	auto start4 = high_resolution_clock::now();
	std::string mode = "betaHighCircle"; //betaHighLune
	bool val = checkInsertSubDsimplex(repsimplex,simpl,distMatrix,beta,root,mode);
	if(val){
			return std::make_pair(true,std::make_pair(simplex,std::make_pair(cc,radius)));			 

	}
	return std::make_pair(false,std::make_pair(simplex,std::make_pair(cc,radius)));			 
}
std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> dwaytreenode:: validatesimplices(dwaytreenode* root,std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> tovalidate,double beta){
   std::set<std::set<std::vector<double>,lexical_compare_points>> validatedSimplices;
   int k;
   std::set<std::vector<double>,lexical_compare_points> accountedpoints;
   std::set<std::vector<double>,lexical_compare_points> totalpoints;
   for(auto x:tovalidate)
	   for(auto y:x.first)
		  totalpoints.insert(y);
   
   while(tovalidate.size()>0){
	    auto simp = *(tovalidate.begin());
	    std::vector<std::vector<double>> simpl(simp.first.begin(), simp.first.end());
		tovalidate.erase(tovalidate.begin());
		auto start4 = high_resolution_clock::now();
		bool val = checkPointInBall(root, simp.second.first,simp.second.second, simpl); 
		auto stop4 = high_resolution_clock::now();
		auto duration4 = duration_cast<microseconds>(stop4 - start4);
		time4.push_back(duration4.count());
	
   		if(!val){
			validatedSimplices.insert(simp.first);
		    std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> remaining;
			for(auto p: tovalidate){
				bool val = true;
				for(auto v:simp.first){
					bool val1 = int(round(utils::vectors_distance(v,p.second.first)*precision))<int(round(p.second.second*precision));
				     if(int(round(utils::vectors_distance(v,p.second.first)*precision))<int(round(p.second.second*precision))){
						val = false;
						break;	
					}
				}
				if(val){
					for(auto v:p.first){
						bool val1 = int(round(utils::vectors_distance(v,simp.second.first)*precision))<int(round(simp.second.second*precision));
						if(int(round(utils::vectors_distance(v,simp.second.first)*precision))<int(round(simp.second.second*precision))){
							val = false;
							break;	
						}
					}
				}
			   if(val){
				remaining.insert(p);
				}
				
			}
		tovalidate.clear();
		tovalidate.insert(remaining.begin(), remaining.end());
		for(auto x:simp.first)
			accountedpoints.insert(x);
		}
   }
   std::vector<std::vector<double>> tp(totalpoints.begin(),totalpoints.end());
	for(auto x : accountedpoints)
		tp.erase(std::remove(tp.begin(), tp.end(), x), tp.end());
	std::set<std::vector<double>,lexical_compare_points> isolatedpoints(tp.begin(),tp.end());
   return std::make_pair(validatedSimplices,isolatedpoints);
}


std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> dwaytreenode::generatePartitionFromIsolatedPoints(dwaytreenode* root,std::set<std::vector<double>,lexical_compare_points> isolatedpoints,int homologydim,double epsilon,double beta){
        std::set<std::set<std::vector<double>,lexical_compare_points>> validones;
        if(isolatedpoints.size()<homologydim+1)
           return std::make_pair(validones,isolatedpoints);
        auto simplices = generateCombinations(isolatedpoints,homologydim);
		std::set<std::vector<double>,lexical_compare_points> accountedpoints;
		std::set<std::vector<double>,lexical_compare_points> remaningpoints;
        for(auto x:simplices){
			double weight = 0;
			for(auto y:x){
				for(auto z:x){
					double dist = utils::vectors_distance(y,z);
					if(dist>weight)
					   weight = dist;
				}
			}
			if(weight <epsilon){
				auto vali = validatesimplex(root,x,homologydim,beta);
				if(vali.first){
				validones.insert(x);
				for(auto z: x)
						accountedpoints.insert(z);
				}
			}
		}
		std::vector<std::vector<double>> tp(isolatedpoints.begin(),isolatedpoints.end());

		for(auto x : accountedpoints)
			tp.erase(std::remove(tp.begin(), tp.end(), x), tp.end());
		std::set<std::vector<double>,lexical_compare_points> isolatedpointsnew(tp.begin(),tp.end());
	
		return std::make_pair(validones,isolatedpointsnew);
}
std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> dwaytreenode:: generateAllSimplicestoCheck1(dwaytreenode* root,std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> propermeshestomerge,std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> newpartition,int homologydim,double beta){
	std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> allmeshes(propermeshestomerge.begin(),propermeshestomerge.end());
	std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> candidateSimplices;
	allmeshes.push_back(newpartition.first);	
    int l;
	for(int k = 0;k<homologydim;k++){
		int otherk = homologydim-k-1;
		int i=0;
		for(auto x:allmeshes){
			auto facetsa = generateCombinationsall(x,k);
			auto face = generateCombinations(newpartition.second,k);
			for(auto v:face)
				facetsa.insert(v);
				
			std::set<std::set<std::vector<double>,lexical_compare_points>> facetsb;
			for(int g=0;g<allmeshes.size();g++){
				if(i!=g){
					auto y = allmeshes[g];
					auto facetsbadd = generateCombinationsall(y,otherk);
					std::set<std::set<std::vector<double>,lexical_compare_points>> simps;
					std::set_union(facetsbadd.begin(), facetsbadd.end(),facetsb.begin(), facetsb.end(),std::inserter(simps, simps.begin()));
					facetsb = simps;
				}
			}
			for(auto x: facetsa){
				for(auto z : facetsb){
					std::set<std::vector<double>,lexical_compare_points> newsimplex;
					std::set_union(x.begin(), x.end(),z.begin(), z.end(),std::inserter(newsimplex, newsimplex.begin()));
					auto simple = validatesimplex(root,newsimplex,homologydim,beta);
					candidateSimplices.insert(simple.second);

				}
			}
			i++;
		}
	}
	
	//Cross Partition simplices for isolated points
	std::set<std::set<std::vector<double>,lexical_compare_points>> candidateSimplices1;
	for(int k = 0;k<homologydim-1;k++){
		int otherk = homologydim-k-2;
		int i=0;
		for(auto x:allmeshes){
			auto facetsa = generateCombinationsall(x,k);	
			std::set<std::set<std::vector<double>,lexical_compare_points>> facetsb;
	
			for(int g=0;g<allmeshes.size();g++){
				if(i!=g){
					auto y = allmeshes[g];
					auto facetsbadd = generateCombinationsall(y,otherk);
					std::set<std::set<std::vector<double>,lexical_compare_points>> simps;
					std::set_union(facetsbadd.begin(), facetsbadd.end(),facetsb.begin(), facetsb.end(),std::inserter(simps, simps.begin()));
					facetsb = simps;
				}
			}
			for(auto x: facetsa){
				for(auto z : facetsb){
					std::set<std::vector<double>,lexical_compare_points> newsimplex;
					std::set_union(x.begin(), x.end(),z.begin(), z.end(),std::inserter(newsimplex, newsimplex.begin()));
					candidateSimplices1.insert(newsimplex);
				}
			}
			i++;
		}
	}
	
	for(auto v:newpartition.second)
	  for(auto tr : candidateSimplices1){
	  		std::set<std::vector<double>,lexical_compare_points> newsimplex(tr.begin(),tr.end());
			newsimplex.insert(v);
			auto simple = validatesimplex(root,newsimplex,homologydim,beta);
			candidateSimplices.insert(simple.second);
		}
			
	return candidateSimplices;
}
std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> dwaytreenode::mergedmesh(dwaytreenode* root, std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> meshestomerge,double beta, int homologydim,double epsilon){
	
	// Coding now only for binary tree, which is suffucient in most practical applications
    std::set<std::vector<double>,lexical_compare_points> isolatedpoints1;
    std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> propermeshestomerge;
    //Collect All the Isolated Points accross all partitions
    for(auto x:meshestomerge){
		for(auto y :x.second){
			isolatedpoints1.insert(y);
		}
		if(x.first.size()>0){
			propermeshestomerge.push_back(x.first);
		}
	}
	//Create a third partition of these Isolated Points by validating any simplex in there
	auto newpartition = generatePartitionFromIsolatedPoints(root,isolatedpoints1,homologydim,epsilon,beta);
	
    if(propermeshestomerge.size()<=0){
		return newpartition;
	}
	 
	auto start2 = high_resolution_clock::now();
	auto propersimplices = generateAllSimplicestoCheck1(root,propermeshestomerge,newpartition,homologydim,beta);
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>(stop2 - start2);
    time2.push_back(duration2.count());
    
	auto start3 = high_resolution_clock::now();
	auto validatedsimplices = validatesimplices(root,propersimplices,beta);	
    auto stop3 = high_resolution_clock::now();
    auto duration3 = duration_cast<microseconds>(stop3 - start3);
    time3.push_back(duration3.count());
	std::vector<std::vector<double>> tp(validatedsimplices.second.begin(),validatedsimplices.second.end());

	
	for(auto x:meshestomerge){
		for(auto y:x.first){
			validatedsimplices.first.insert(y);
			for(auto t:y)
				tp.erase(std::remove(tp.begin(), tp.end(), t), tp.end());
		}
	}
	std::set<std::vector<double>,lexical_compare_points> isolatedpointsnew(tp.begin(),tp.end());
	validatedsimplices.second = isolatedpointsnew;
	for(auto x:newpartition.first)
		validatedsimplices.first.insert(x);
		
    return validatedsimplices;
	/*
	auto propersimplices = generateAllSimplicestoCheck(meshestomerge,homologydim);
    auto validatedsimplices = validatesimplices(root,propersimplices,beta);	
    for(auto x:meshestomerge[0].first){
		validatedsimplices.insert(x);
	}
	for(auto x:meshestomerge[1].first){
		validatedsimplices.insert(x);
	}
	if(isolatedpoints1.size()<=homologydim){
	    return std::make_pair(validatedsimplices,isolatedpoints1);
	}
	else{
		auto isolatedmesh = generateCombinations(isolatedpoints1,homologydim);
		std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> abc;
		std::set<std::vector<double>,lexical_compare_points> iso;
		abc.push_back(std::make_pair(validatedsimplices,iso));
		abc.push_back(std::make_pair(isolatedmesh,iso));
		auto propersimplices1 = generateAllSimplicestoCheck(abc,homologydim);
		std::set<unsigned> repsimplex;
		for(int i=0;i<homologydim+1;i++)
				repsimplex.insert(i);
				
		for(auto tp : isolatedmesh){
			    int i=0,j=0;
			    std::vector<std::vector<double>> simp(tp.begin(), tp.end());
				std::vector<std::vector<double>> distMatrix(simp.size(),std::vector<double>(simp.size()));
		
				for(auto y:simp){
					j=0;
					for(auto z:simp){
						distMatrix[i][j] = utils::vectors_distance(y,z);
						j++;
					}
					i++;
				} 
		
			auto cc =  utils::circumCenter(repsimplex,simp);
			auto radius = sqrt(utils::circumRadius(repsimplex,&distMatrix));
		 
		    propersimplices1.insert(std::make_pair(tp,std::make_pair(cc,radius)));
		}
		auto validatedsimplices1 = validatesimplices(root,propersimplices1,beta);
		for(auto x:validatedsimplices1)
			validatedsimplices.insert(x);
		for(auto g: validatedsimplices1	){
			for(auto d:g){
				isolatedpoints1.erase(d);
			}
		}

	   return std::make_pair(validatedsimplices,isolatedpoints1);
	}

  return std::make_pair(validatedsimplices,isolatedpoints1);
*/
	    
	/*
	std::vector<std::vector<double>> isolatedpoints(isolatedpoints1.begin(), isolatedpoints1.end());
    std::set<std::set<std::vector<double>,lexical_compare_points>> simplicestocheck1;
	std::set<std::set<std::vector<double>,lexical_compare_points>> simplicestocheck;
	std::set<std::vector<double>,lexical_compare_points> newisolates;
	if(isolatedpoints.size()<=homologydim){
		newisolates = isolatedpoints1;
	}
	else{
		simplicestocheck = generateCombinations(isolatedpoints1,homologydim);
		simplicestocheck1 = simplicestocheck;
	}
	int i=0;
	std::set<std::set<std::vector<double>,lexical_compare_points>> validlist;
	for(auto x:meshestomerge){
		std::set<std::vector<double>,lexical_compare_points> pointsinotherpartition;
		for(int g=0;g<meshestomerge.size();g++){
			if(i!=g){
				auto y = meshestomerge[g];
				auto newpts = generatePoints(y);
				std::set<std::vector<double>,lexical_compare_points> pnts;
				std::set_union(newpts.begin(), newpts.end(),pointsinotherpartition.begin(), pointsinotherpartition.end(),std::inserter(pnts, pnts.begin()));
				pointsinotherpartition = pnts;
			}
		}
		auto newsimplices = generateNewSimplices(x.first,pointsinotherpartition);
		for(auto t:x.first)
			validlist.insert(t);
		std::set<std::set<std::vector<double>,lexical_compare_points>> validlistmerged;
		std::set_union(newsimplices.begin(), newsimplices.end(),simplicestocheck.begin(), simplicestocheck.end(),std::inserter(validlistmerged, validlistmerged.begin()));
		simplicestocheck = validlistmerged;
	
		i++;
	}
	std::set<std::vector<double>,lexical_compare_points> pointsinotherpartition1;
	for(int g=0;g<meshestomerge.size();g++){
				auto y = meshestomerge[g];
				auto newpts1 = generatePoints(y);		
				std::set<std::vector<double>,lexical_compare_points> pntstemp;
				std::set_union(newpts1.begin(), newpts1.end(),pointsinotherpartition1.begin(), pointsinotherpartition1.end(),std::inserter(pntstemp, pntstemp.begin()));
				pointsinotherpartition1 = pntstemp;

	}
    std::vector<std::vector<double>> tp(pointsinotherpartition1.begin(),pointsinotherpartition1.end());

	for(auto x : isolatedpoints1)
		tp.erase(std::remove(tp.begin(), tp.end(), x), tp.end());
	std::set<std::vector<double>,lexical_compare_points> pintpart(tp.begin(),tp.end());
	auto newsimplices1 = generateNewSimplices(simplicestocheck1,pintpart);
	
	
	std::set<std::set<std::vector<double>,lexical_compare_points>> validlistmerged1;
	std::set_union(newsimplices1.begin(), newsimplices1.end(),simplicestocheck.begin(), simplicestocheck.end(),std::inserter(validlistmerged1, validlistmerged1.begin()));
	simplicestocheck = validlistmerged1;
	* 
	auto newvalidlist = filterValidSimplices(root,simplicestocheck,beta);
	std::set<std::set<std::vector<double>,lexical_compare_points>> finalvalidlist = validlist;
	std::set<std::set<std::vector<double>,lexical_compare_points>> fvl;
	std::set_union(finalvalidlist.begin(), finalvalidlist.end(),newvalidlist.first.begin(), newvalidlist.first.end(),std::inserter(fvl, fvl.begin()));

    std::set<std::vector<double>,lexical_compare_points> newisolateupdated;
	std::set_union(newisolates.begin(), newisolates.end(),newvalidlist.second.begin(), newvalidlist.second.end(),std::inserter(newisolateupdated, newisolateupdated.begin()));
	
	std::vector<std::vector<double>> kp(isolatedpoints.begin(),isolatedpoints.end());

	for(auto x : newisolateupdated)
		kp.erase(std::remove(kp.begin(), kp.end(), x), kp.end());
	
	std::vector<std::vector<double>> kp1(isolatedpoints.begin(),isolatedpoints.end());

	for(auto x : kp)
		kp1.erase(std::remove(kp1.begin(), kp1.end(), x), kp1.end());
	std::set<std::vector<double>,lexical_compare_points> isolatedpo(kp1.begin(),kp1.end());
	

	return std::make_pair(fvl,newisolateupdated);
	*/
}
std::vector<std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>>> dwaytreenode::computeSpaceFabricNeighbourhood(dwaytreenode* root,std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> properMeshes,double epsilon,int homologydim){
 	
    std::vector<std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>>> filteredSimplices;
    int i=0;
    for(auto mesh: properMeshes){
		int j=0;
		double d = 0;
		int vec1;
		 if(i>=root->directionVectors.size())
			vec1 = 0;
	    else 
 	        vec1 = i;
 		for(auto x:root->coordinates)
			d -=x*root->directionVectors[vec1][j++];
    
		std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> allSimplices(homologydim+1);
		for(auto properMesh : mesh){
			int k=0;
			std::set<std::vector<double>,lexical_compare_points> filteredpts;
			for(auto x : properMesh){
				    int vec;
				    if(i>=root->directionVectors.size())
				       vec = 0;
				    else 
				       vec = i;
					double distancefromSplittingPlane = utils::distanceFromHyperplane(x,root->directionVectors[vec], d);
					if(distancefromSplittingPlane<epsilon){
						k++;
						filteredpts.insert(x);

					}
			}
			if(k>0)	
				allSimplices[k-1].insert(filteredpts);
		}
		filteredSimplices.push_back(allSimplices);
		i++;	
	}	
	
    return filteredSimplices;
}

std::set<std::set<std::vector<double>,lexical_compare_points>> dwaytreenode:: generateKSimplices(std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> stichNeighboorhood,int k,int homologydim){
    std::set<std::set<std::vector<double>,lexical_compare_points>> Ksimplices(stichNeighboorhood[k].begin(),stichNeighboorhood[k].end());
    
    for(int i=k+1;i<homologydim+1;i++){
	     auto ksim = generateCombinationsall(stichNeighboorhood[i],k);
	     for(auto g : ksim)
	          Ksimplices.insert(g);
	}
  
	return Ksimplices;
}


std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> dwaytreenode:: generateSimplicesToConsider(dwaytreenode* root,std::vector<std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>>> stichNeighboorhood,std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> propermeshestomerge,std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> newpartition,int homologydim,double epsilon,double beta){
	std::set<std::pair<std::set<std::vector<double>,lexical_compare_points>,std::pair<std::vector<double>,double>>,comp_by_radius> simpliceToConsider;
	for(int k = 0;k<homologydim;k++){
		int otherk = homologydim-k-1;
		for(int i = 0;i<stichNeighboorhood.size();i++){
			auto ksimplices = generateKSimplices(stichNeighboorhood[i],k,homologydim); 
			auto face = generateCombinations(newpartition.second,k);
			for(auto v:face)
				ksimplices.insert(v);	
            for(auto ksim : ksimplices){
		
				for(int g = 0;g<stichNeighboorhood.size();g++){
					if(i!=g){
						std::set<std::vector<double>> intersectionregion;
						bool first = true;
						for(auto pt :ksim){
							if(first)
								intersectionregion = pointInBall(root,pt,epsilon);
							else{
								auto intersectionregionnext =  pointInBall(root,pt,epsilon); //find regionofinterestinotherpartion
								for(auto d: intersectionregionnext)
									if(intersectionregion.find(d)==intersectionregion.end())
										intersectionregion.erase(d);
							}
						      first = false;
						}
						auto otherksimplices = generateKSimplices(stichNeighboorhood[g],otherk,homologydim); 
						for(auto otherksim : otherksimplices){
							int count =0;
							for(auto othsim : otherksim){
								if(intersectionregion.find(othsim) !=intersectionregion.end())
								   count++;
								}
							if(count==otherksim.size()){
								// Mergetobecheckedforvalidsimplices
								std::set<std::vector<double>,lexical_compare_points> newsimplex;
								std::set_union(ksim.begin(), ksim.end(),otherksim.begin(), otherksim.end(),std::inserter(newsimplex, newsimplex.begin()));
								if(newsimplex.size()==homologydim+1){
									auto simple = validatesimplex(root,newsimplex,homologydim,beta);
									simpliceToConsider.insert(simple.second);
								}

							}
						}
					}
				}
			}
		}
	}
	std::set<std::set<std::vector<double>,lexical_compare_points>> simpliceToConsider1;
	for(int k = 0;k<homologydim-1;k++){
		int otherk = homologydim-k-2;
		for(int i = 0;i<stichNeighboorhood.size();i++){
			auto ksimplices = generateKSimplices(stichNeighboorhood[i],k,homologydim); 
	        for(auto ksim : ksimplices){
				for(int g = 0;g<stichNeighboorhood.size();g++){
					if(i!=g){
						std::set<std::vector<double>> intersectionregion;
						bool first = true;
						for(auto pt :ksim){
							if(first)
								intersectionregion = pointInBall(root,pt,epsilon);
							else{
								auto intersectionregionnext =  pointInBall(root,pt,epsilon); //find regionofinterestinotherpartion
								for(auto d: intersectionregionnext)
									if(intersectionregion.find(d)==intersectionregion.end())
										intersectionregion.erase(d);
							}
						      first = false;
						}
						auto otherksimplices = generateKSimplices(stichNeighboorhood[g],otherk,homologydim); 
						for(auto otherksim : otherksimplices){
							int count =0;
							for(auto othsim : otherksim){
								if(intersectionregion.find(othsim) !=intersectionregion.end())
								   count++;
								}
							if(count==otherksim.size()){
								// Mergetobecheckedforvalidsimplices
								std::set<std::vector<double>,lexical_compare_points> newsimplex;
								std::set_union(ksim.begin(), ksim.end(),otherksim.begin(), otherksim.end(),std::inserter(newsimplex, newsimplex.begin()));
								if(newsimplex.size()==homologydim){			
									simpliceToConsider1.insert(newsimplex);
								}
							}
						}
					}
				}
			}
		}
	}
	
	if(stichNeighboorhood.size()==1){
		auto ksimplices = generateKSimplices(stichNeighboorhood[0],homologydim-1,homologydim);
		for(auto pt :newpartition.second){
			for(auto subsim :ksimplices){
				std::set<std::vector<double>,lexical_compare_points> newsimplex(subsim.begin(),subsim.end());
				newsimplex.insert(pt);
				if(newsimplex.size()==homologydim+1){				
				auto simple = validatesimplex(root,newsimplex,homologydim,beta);
				simpliceToConsider.insert(simple.second);
				}
			}
		}
	}	
	
		
	for(auto v:newpartition.second)
	  for(auto tr : simpliceToConsider1){
	  		std::set<std::vector<double>,lexical_compare_points> newsimplex(tr.begin(),tr.end());
			newsimplex.insert(v);
				if(newsimplex.size()==homologydim+1){			
					auto simple = validatesimplex(root,newsimplex,homologydim,beta);
					simpliceToConsider.insert(simple.second);
				}
		}
		
	return simpliceToConsider;
}

std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> dwaytreenode::stichMesh(dwaytreenode* mainroot,dwaytreenode* root, std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> meshestomerge,double beta, int homologydim,double epsilon){
	
	
	// Coding now only for binary tree, which is suffucient in most practical applications
    std::set<std::vector<double>,lexical_compare_points> isolatedpoints1;
    std::set<std::vector<double>,lexical_compare_points> isolatedpointsnotinrange;
    std::set<std::vector<double>,lexical_compare_points> isolatedpointsinrange;

    std::vector<std::set<std::set<std::vector<double>,lexical_compare_points>>> propermeshestomerge;
    
    //Collect All the Isolated Points accross all partitions
    int i=0;
    for(auto x:meshestomerge){
		int j=0;
		double d = 0;
		for(auto x:root->coordinates)
			d -=x*root->directionVectors[i][j++];
		for(auto y :x.second){
			isolatedpoints1.insert(y);
			if(utils::distanceFromHyperplane(y,root->directionVectors[i], d)<epsilon)
				isolatedpointsinrange.insert(y);
			else
				isolatedpointsnotinrange.insert(y);
		}
		if(x.first.size()>0){
			propermeshestomerge.push_back(x.first);
		}
		i++;
	}

	//Create a third partition of these Isolated Points by validating any simplex in there
	
	auto newpartition = generatePartitionFromIsolatedPoints(mainroot,isolatedpointsinrange,homologydim,epsilon,beta);
	for(auto x:isolatedpointsnotinrange){
		newpartition.second.insert(x);
	}
    if(propermeshestomerge.size()<=0){
		return newpartition;
	}

	if(newpartition.first.size()>0)
		propermeshestomerge.push_back(newpartition.first);

	auto stichNeighboorhood = computeSpaceFabricNeighbourhood(root,propermeshestomerge,epsilon,homologydim);	
			
    auto properSimplices = generateSimplicesToConsider(root,stichNeighboorhood,propermeshestomerge,newpartition,homologydim,epsilon,beta);
    
    auto validatedsimplices = validatesimplices(mainroot,properSimplices,beta);	

	std::vector<std::vector<double>> tp(validatedsimplices.second.begin(),validatedsimplices.second.end());

	for(auto x:meshestomerge){
		for(auto y:x.first){
			validatedsimplices.first.insert(y);
			for(auto t:y)
				tp.erase(std::remove(tp.begin(), tp.end(), t), tp.end());
		}
	}
	std::set<std::vector<double>,lexical_compare_points> isolatedpointsnew(tp.begin(),tp.end());
	validatedsimplices.second = isolatedpointsnew;
	for(auto x:newpartition.first)
		validatedsimplices.first.insert(x);
		
    return validatedsimplices;
}

std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> dwaytreenode::meshGeneration(dwaytreenode* mainroot,dwaytreenode* root, double beta, int homologydim, double epsilon){
	std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>> mesh;		
    if (root == nullptr){
        return mesh;
	}
    if(root->children.size()<=0){
		mesh.second.insert(root->coordinates);
		return mesh;
	}
 	std::vector<std::pair<std::set<std::set<std::vector<double>,lexical_compare_points>>,std::set<std::vector<double>,lexical_compare_points>>> meshestomerge;		
    
    
    for(auto child:root->children){
		meshestomerge.push_back(meshGeneration(mainroot,child,beta,homologydim,epsilon));
	}
	return stichMesh(mainroot,root,meshestomerge,beta,homologydim,epsilon);

//	return mergedmesh(mainroot,meshestomerge,beta,homologydim,epsilon);

}

dwaytreenode* dwaytreenode::buildDwayTree(std::vector<std::vector<double>> data,int direction,std::string mode,int directionCount){
	std::vector<double> centroid;
	std::vector<std::vector<double>> coords = utils::transpose(data);
	for(auto coordDim : coords){
		centroid.push_back(utils::getAverage(coordDim));
	}
	double radius = 0;
	for(auto x: data){
		double distance = utils::vectors_distance(x,centroid);
		if(distance>radius)
		    radius = distance;
	}
	if(data.size()<=1){
		return new dwaytreenode(centroid,radius,direction);
	}
    //**********************Partition the data into d+1 buckets*******************
    // Project all the points on unut hypersphere and run kmeans clustering to get required directions vectors and partiotions.
    // Will use existing kmeans++ code to get it working here.
    std::vector<unsigned> partitions;
    std::vector<std::vector<double>> dv = referenceHypertetrhedron;
    std::pair<std::vector<std::vector<double>>,std::vector<unsigned>> parts;
    if(mode == "Default" || mode == "default"){
		for(auto A : data){
			std::vector<double> A_vec;
			int k =0;
			for(auto a:A){
				A_vec.push_back(a-centroid[k]);
				k++;
			}
			int direction = -1;
			int assignedpartion = 0;
			int maxvalue = 0;
			for(auto B : referenceHypertetrhedron){
				direction = direction+1;
				double cosine =  utils::cosine_similarity(A_vec,B)*(180/M_PI);
				if(maxvalue<cosine){
					maxvalue = cosine;
					assignedpartion = direction;
				}
			}
			partitions.push_back(assignedpartion);
		}
	}
	else if(mode == "nary" || mode == "Nary"){
		parts = utils::computeDirectionVectors(data,centroid,directionCount);
		partitions = parts.second;
		dv = parts.first;
	}
	
	dwaytreenode* root = new dwaytreenode(centroid,radius,direction,dv);

	std::vector<std::vector<std::vector<double>>> dwaypartition(directionCount, std::vector<std::vector<double>>(0, std::vector<double>(0)));
	int pts = 0;
	for(auto x : partitions){
		dwaypartition[x].push_back(data[pts]);
		pts = pts+1;
	}		
	//*****************************Partioning Done********************************
	direction =0;
	// Recursively build the tree for each partiotion. Stop recursion when the split size is 1 or 0.
	for(auto part : dwaypartition){
		if(part.size()>1){
			dwaytreenode *child = root->buildDwayTree(part,direction,mode,directionCount);
			child->parent = root;
			root->children.push_back(child);
		}
		else if(part.size()==1){
			dwaytreenode *child = new dwaytreenode(part[0],0,direction);
			child->parent = root;
			root->children.push_back(child);
		}
		direction++;
	}	
    return root;
    
}
 
void dwaytreenode::pathToCell(std::vector<dwaytreenode*> &path,dwaytreenode* croot, dwaytreenode* node){
	int direction = -1;
	int assignedpartition = 0;
	int maxvalue = 0;
	path.push_back(croot);
	if(croot==node)
		return;
	for(auto B : referenceHypertetrhedron){
		direction = direction+1;
		std::vector<double> A_vec;
		int k =0;
		for(auto a:croot->coordinates){
			A_vec.push_back(node->coordinates[k]-a);
			k++;
		}
		double cosine =  utils::cosine_similarity(A_vec,B)*(180/M_PI);
		if(maxvalue<cosine){
			maxvalue = cosine;
			assignedpartition = direction;
		}
	}
   for(auto child : croot->children){
	   if(child->parenttochilddirection==assignedpartition){
			pathToCell(path,child,node);
		}
   }
   return;
   
}

std::vector<std::vector<double>>  dwaytreenode::findCellBoundingPolytope(dwaytreenode* root, dwaytreenode* node){
      std::vector<dwaytreenode*> path;
      pathToCell(path,root,node);
		
	  std::vector<std::vector<double>> boundingVertices(dim+1, std::vector<double>(0));
	  
	  for(auto x:path){
		  if(x->parenttochilddirection !=-1){
				boundingVertices[x->parenttochilddirection] = x->parent->coordinates;
		  }
	  }  
	  return boundingVertices;	
}


void dwaytreenode::printTree(dwaytreenode* root){
	 dwaytreenode* temp = root;
     if (temp == nullptr) {
        std::cout << "Tree empty" << std::endl;
        return;
    }
	for(auto child : temp->children)
		printTree(child);
	
	for(auto x : temp->coordinates)
     std::cout<<x<<" ";
    std::cout<<"\n";
    
}

void makeCombiUtil(std::vector<std::vector<int> >& ans,std::vector<int>& tmp, int n, int left, int k)
{
    // Pushing this vector to a vector of vector
    if (k == -1) {
        ans.push_back(tmp);
        return;
    }
 
    // i iterates from left to n. First time
    // left will be 1
    for (int i = left; i <= n; ++i)
    {
        tmp.push_back(i);
        makeCombiUtil(ans, tmp, n, i + 1, k - 1);
 
        // Popping out last inserted element
        // from the vector
        tmp.pop_back();
    }
}
 
// Prints all combinations of size k of numbers
// from 1 to n.
std::vector<std::vector<int> > makeCombi(int n, int k)
{
    std::vector<std::vector<int> > ans;
    std::vector<int> tmp;
    makeCombiUtil(ans, tmp, n-1, 0, k-1);
    return ans;
}

dwaytreenode* dwaytreenode:: findNearestNeighbor(dwaytreenode* root, std::vector<double> pt){
	int direction = -1;
	int assignedpartion = 0;
	int maxvalue = 0;
	int k =0;

	if(root->children.size()<=1)
		return root;
	std::vector<double> A_vec;
		for(auto a:pt){
			A_vec.push_back(a-root->coordinates[k]);
			k++;
		}
	for(auto B : root->directionVectors){
		direction = direction+1;
		double cosine =  utils::cosine_similarity(A_vec,B)*(180/M_PI);
		if(maxvalue<cosine){
			maxvalue = cosine;
			assignedpartion = direction;
		}
	}
	dwaytreenode* temp = root;
	for(auto x:root->children){
	if(x!=nullptr)
		if(x->parenttochilddirection == assignedpartion)
			temp = findNearestNeighbor(x,pt);	
	}

	dwaytreenode* best = temp;
	double bestradius = utils::vectors_distance(pt,temp->coordinates);
	
	for(auto x:root->children){
		double distance = utils::vectors_distance(pt,x->coordinates);
		if(distance<bestradius+x->radius){
			if(x->parenttochilddirection != assignedpartion){
				temp = findNearestNeighbor(x,pt);
				double radius = utils::vectors_distance(pt,temp->coordinates);
				if(bestradius > radius){
				    bestradius = radius;
				    best = temp;
				}
			}
		}
	}
	
	return best;
	
}
	
bool dwaytreenode::checkPointInBall(dwaytreenode* root, std::vector<double> pt,double radius, std::vector<std::vector<double>> omissions)
{   
	if(root->children.size()==0)
		if(utils::vectors_distance(pt,root->coordinates) < radius)
			return std::find(omissions.begin(),omissions.end(),root->coordinates)==omissions.end() ? true : false;
	for(auto child_node : root->children)
		if(utils::vectors_distance(pt,child_node->coordinates)<(radius+child_node->radius))
			if(checkPointInBall(child_node, pt, radius,omissions))
				return true;
	return false;
}

std::set<std::vector<double>> dwaytreenode::pointInBall(dwaytreenode* root, std::vector<double> pt,double radius)
{
	std::set<std::vector<double>> point_in_Ball;
	if(root->children.size()==0)
	{
		if(utils::vectors_distance(pt,root->coordinates) < radius)
			point_in_Ball.insert(root->coordinates);
		return point_in_Ball;
	}
	for(auto child_node : root->children)
		if(utils::vectors_distance(pt,child_node->coordinates)<(radius+child_node->radius))
			for(auto point:pointInBall(child_node, pt, radius))
				point_in_Ball.insert(point);
	return point_in_Ball;
}


bool dwaytreenode:: checkInsertSubDsimplex(std::set<unsigned> dsimplex,std::vector<std::vector<double>> data,std::vector<std::vector<double>> distMatrix,double beta,dwaytreenode* tree,std::string mode){
	std::vector<std::vector<double>> neighborsfinalLune;
	std::vector<std::vector<double>> neighborsfinalCircle;
        
	bool intersectionCircle= false;
	bool intersectionLune = false;
	bool unionCircle= false;

	if(beta <0)
		exit(0);
	else if(beta==0)
		return true;
    else if(beta <1)
		intersectionCircle = true;
	else if(beta >1){
		if(mode == "betaHighLune")
		     intersectionLune = true;
		else
			 unionCircle = true;
	}
	else if(beta==1){
		intersectionLune=true;
		intersectionCircle=true;
		unionCircle = true;
	}
	
	if(mode == "betaHighCircle"){
		if(beta<1)
			beta=1/beta;
		std::set<unsigned> simplex(dsimplex.begin(),dsimplex.end());
	    std::vector<double> circumCenter;
	   	if(simplex.size()>2)
			circumCenter = utils::circumCenter(simplex,data);
		else if(simplex.size()==2){
 			auto first = simplex.begin();
			std::vector<double> R;
			std::vector<double> A = data[*first];
			std::advance(first, 1);
			std::vector<double> B = data[*first];
			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
			circumCenter = R;
   		}
                
		double circumRadius;
        circumRadius = pow(utils::vectors_distance(circumCenter, data[0]),2);
		std::vector<std::vector<double>> neighbors;
		std::vector<std::vector<std::vector<double>>> neighborsCircleIntersection;   
		std::vector<std::vector<double>> mat;
   
		for(auto x : simplex)
			mat.push_back(data[x]);
		mat.push_back(circumCenter);
		
		auto hpcoff = utils::nullSpaceOfMatrix(simplex,data,circumCenter,sqrt(circumRadius),true);

        std::vector<std::vector<double>> refbetaCenters ;
		refbetaCenters = utils::betaCentersCalculation(hpcoff.first, beta, sqrt(circumRadius),circumCenter);
		
		refbetaCenters = utils::computePCAInverse(mat,refbetaCenters,hpcoff.second);
		double betaRadius = utils::vectors_distance(refbetaCenters[0], data[0]);
        std::set<std::vector<double>> neighbors1 = pointInBall(tree,refbetaCenters[0], betaRadius); //All neighbors in epsilon-ball
		for(auto t :data)
			neighbors1.erase(t);
		std::set<std::vector<double>> neighbors2 = pointInBall(tree,refbetaCenters[1], betaRadius); //All neighbors in epsilon-ball
		if(intersectionCircle==true){
			std::vector<std::vector<double>> v1(std::min(neighbors1.size(),neighbors2.size()));
			std::vector<std::vector<double>>::iterator it1;
			it1=std::set_intersection (neighbors1.begin(), neighbors1.end(), neighbors2.begin(), neighbors2.end(), v1.begin());
			v1.resize(it1-v1.begin()); 
			std::vector<std::vector<double>> k1(v1.begin(),v1.end());
			neighborsfinalCircle = k1;
		}
		else if(unionCircle ==true){
			std::vector<std::vector<double>> v1(neighbors1.size()+neighbors2.size());
			std::vector<std::vector<double>>::iterator it1;				
			it1=std::set_union(neighbors1.begin(), neighbors1.end(), neighbors2.begin(), neighbors2.end(), v1.begin());
			v1.resize(it1-v1.begin());
			std::vector<std::vector<double>> k1(v1.begin(),v1.end()); 
			neighborsfinalCircle = k1;
		}
	}

	else if(mode == "betaHighLune"){
		std::set<unsigned> simplex(dsimplex.begin(),dsimplex.end());
       		std::vector<double> circumCenter;
		std::vector<double> circumCenterfaces;
		std::vector<double> circumCenterfaces1;
		if(simplex.size()>2)
			circumCenter = utils::circumCenter(simplex,data);
		else if(simplex.size()==2){
			auto first = simplex.begin();
			std::vector<double> R;
			std::vector<double> A = data[*first];
			std::advance(first, 1);
			std::vector<double> B = data[*first];
			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
			circumCenter = R;
		}
		double circumRadius;
        circumRadius = pow(utils::vectors_distance(circumCenter, data[0]),2);

		bool first = true;
		for (auto x : simplex){
			std::vector<double> betaCenter;
			for(unsigned y =0 ;y< data[0].size();y++)
				betaCenter.push_back(beta*circumCenter[y] + (1-beta)*data[x][y]);
			double betaRadius = beta*sqrt(circumRadius);     
		
			std::set<std::vector<double>> neighbors1faces1 = pointInBall(tree,betaCenter, betaRadius); //All neighbors in epsilon-ball
			neighbors1faces1.erase(data[x]);
			if(!first){
				if(intersectionLune==true){
					std::vector<std::vector<double>> v1(std::min(neighborsfinalLune.size(),neighbors1faces1.size()));
					std::vector<std::vector<double>>::iterator it1;
					it1=std::set_intersection (neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1-v1.begin()); 
					std::vector<std::vector<double>> k1(v1.begin(),v1.end());
					neighborsfinalLune = k1;
				}
				else{
					std::vector<std::vector<double>> v1(neighborsfinalLune.size()+neighbors1faces1.size());
					std::vector<std::vector<double>>::iterator it1;				
					it1=std::set_union(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1-v1.begin()); 
					std::vector<std::vector<double>> k1(v1.begin(),v1.end());
					neighborsfinalLune=k1;
				}
			}	
			else{
				std::vector<std::vector<double>> k1(neighbors1faces1.begin(),neighbors1faces1.end());
				neighborsfinalLune = k1;
			}
			first= false;
		}
	}
	if(mode == "betaHighLune" && neighborsfinalLune.size() == 0)
		return true;
	if(mode == "betaHighCircle" && neighborsfinalCircle.size() == 0){
		return true;
		
	}
return false;
}

