
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

#define MAXRADIUS 9999

std::vector<std::vector<double>> referenceHypertetrhedron;
int dim;

void makeCombiUtil(std::vector<std::vector<int> >& ans,std::vector<int>& tmp, int n, int left, int k);
std::vector<std::vector<int> > makeCombi(int n, int k);


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
	std::vector<std::vector<double>> pointInBall(dwaytreenode* root, std::vector<double>,double);
	std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> meshGeneration(dwaytreenode* mainroot,dwaytreenode* root, double beta, int homologydim);
	std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> mergedmesh(dwaytreenode* root, std::vector<std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>>> meshstomerge, double beta, int homologydim);
	std::pair<std::vector<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> filterValidSimplices(dwaytreenode* root,std::set<std::vector<std::vector<double>>> simplicestocheck,double beta);
	std::set<std::vector<std::vector<double>>> generateNewSimplices(std::set<std::vector<std::vector<double>>> simplices, std::vector<std::vector<double>> points);
	std::set<std::vector<double>> generatePoints(std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> partition);
	std::set<std::vector<std::vector<double>>> generateCombinations(std::vector<std::vector<double>> isolatedpoints,double homologydim);
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


std::pair<std::vector<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> dwaytreenode::filterValidSimplices(dwaytreenode* root,std::set<std::vector<std::vector<double>>> simplicestocheck,double beta){
	std::vector<std::vector<std::vector<double>>> validlist;
	std::set<unsigned> simplex;
	if(simplicestocheck.size()>0)
	
		for(int i=0;i<(*simplicestocheck.begin()).size();i++)
			simplex.insert(i);
	std::set<std::vector<double>> pointsAccounted;
	std::set<std::vector<double>> totalPoints;
	for(auto x:simplicestocheck){
		for(auto y:x){
			for(auto z:y){
				std::cout<<z<<",";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
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
  		auto cc =  utils::circumCenter(simplex,x);
		auto radius = sqrt(utils::circumRadius(simplex,&distMatrix));
		std::cout<<"\n"<<" "<<radius<<"\n";
		for(auto g:cc)
			std::cout<<g<<" ";
		std::cout<<"\n";
		if(!checkPointInBall(root, cc,radius, x)){
			validlist.push_back(x);
			for(auto pt:x){
				pointsAccounted.insert(pt);
			}
			std::cout<<"Valid\n";
		}	
		for(auto pt:x){
				totalPoints.insert(pt);
			}
			
	}
	std::vector<std::vector<double>> isolatedpoints;
	set_difference(begin(totalPoints), end(totalPoints),  begin(pointsAccounted), end(pointsAccounted), inserter(isolatedpoints, end(isolatedpoints))); 
	
	//Need to filter out list based on beta inclusion rule.
	//Need to report isolated points also
	return std::make_pair(validlist,isolatedpoints);
}

std::set<std::vector<std::vector<double>>> dwaytreenode:: generateNewSimplices(std::set<std::vector<std::vector<double>>> simplices,std::vector<std::vector<double>> points){
	std::set<std::vector<std::vector<double>>> newsimplices;

	for(auto x: simplices){
		auto facets = generateCombinations(x,x.size()-2);
		for(auto z : facets){
			for(auto y : points){
				z.push_back(y);
				newsimplices.insert(z);
				z.pop_back();
			}
		}
	}
	//It will generate simplices by taking facets from A with each points in B
	return newsimplices;
}

std::set<std::vector<double>> dwaytreenode:: generatePoints(std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> partition){
	std::set<std::vector<double>> points;
	
	for(auto x: partition.second)
		points.insert(x);
	
	for(auto y: partition.first)
		for(auto z:y)
			points.insert(z);
	
	//It will degenerate partition to its points
	return points;
}

std::set<std::vector<std::vector<double>>> dwaytreenode:: generateCombinations(std::vector<std::vector<double>> isolatedpoints,double homologydim){
    std::set<std::vector<std::vector<double>>> simplices;
    auto simp = makeCombi(isolatedpoints.size(), homologydim+1);
    for(auto x : simp){
		std::vector<std::vector<double>> sim;
		for(auto y:x){
			sim.push_back(isolatedpoints[y]);
		}
		simplices.insert(sim);
	}
    //generate simplices from isolate points.
    return simplices;
}

std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> dwaytreenode::mergedmesh(dwaytreenode* root, std::vector<std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>>> meshestomerge,double beta, int homologydim){
	
	// Coding now only for binary tree, which is suffucient in most practical applications
    std::vector<std::vector<double>> isolatedpoints;
    for(auto x:meshestomerge){
		for(auto y :x.second){
			isolatedpoints.push_back(y);
		}
	}
	std::set<std::vector<std::vector<double>>> simplicestocheck;
	std::vector<std::vector<double>> newisolates;
	
	if(isolatedpoints.size()<=homologydim){
		newisolates = isolatedpoints;
	}
	else{
		simplicestocheck = generateCombinations(isolatedpoints,homologydim);
	}
	int i=0;
	std::set<std::vector<std::vector<double>>> validlist;
	for(auto x:meshestomerge){
		std::vector<std::vector<double>> pointsinotherpartition;
		for(int g=0;g<meshestomerge.size();g++){
			if(i!=g){
				auto y = meshestomerge[g];
				auto newpts = generatePoints(y);		
				pointsinotherpartition.insert(pointsinotherpartition.end(), newpts.begin(), newpts.end());
			}
		}
		auto newsimplices = generateNewSimplices(x.first,pointsinotherpartition);
		for(auto t:x.first)
			validlist.insert(t);
		std::set<std::vector<std::vector<double>>> validlistmerged;
		std::set_union(newsimplices.begin(), newsimplices.end(),simplicestocheck.begin(), simplicestocheck.end(),std::inserter(validlistmerged, validlistmerged.begin()));
		simplicestocheck = validlistmerged;
		i++;
	}
	auto newvalidlist = filterValidSimplices(root,simplicestocheck,beta);
	std::set<std::vector<std::vector<double>>> finalvalidlist = validlist;
	std::set<std::vector<std::vector<double>>> fvl;
	std::set_union(finalvalidlist.begin(), finalvalidlist.end(),newvalidlist.first.begin(), newvalidlist.first.end(),std::inserter(fvl, fvl.begin()));
    newisolates.insert(newisolates.end(),newvalidlist.second.begin(),newvalidlist.second.end());
	return std::make_pair(fvl,newisolates);
}


std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> dwaytreenode::meshGeneration(dwaytreenode* mainroot,dwaytreenode* root, double beta, int homologydim){
	std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>> mesh;		
    if (root == nullptr)
        return mesh;
    if(root->children.size()<=0){
		mesh.second.push_back(root->coordinates);
		return mesh;
	}
 	std::vector<std::pair<std::set<std::vector<std::vector<double>>>,std::vector<std::vector<double>>>> meshestomerge;		
    
    
    for(auto child:root->children){
		meshestomerge.push_back(meshGeneration(root,child,beta,homologydim));
	}
	for(auto mesh:meshestomerge){
	std::cout<<"simplices::\n";
    for(auto x:mesh.first){
		for(auto y:x){
			for(auto z:y){
				std::cout<<z<<",";
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
	std::cout<<"Points\n";
	for(auto x:mesh.second){
		for(auto y:x){
			std::cout<<y<<",";
		}
		std::cout<<std::endl;
	}    
	}
	
	int k;
	std::cin>>k;
	return mergedmesh(mainroot,meshestomerge,beta,homologydim);

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

std::vector<std::vector<double>> dwaytreenode::pointInBall(dwaytreenode* root, std::vector<double> pt,double radius)
{
	std::vector<std::vector<double>> point_in_Ball;
	if(root->children.size()==0)
	{
		if(utils::vectors_distance(pt,root->coordinates) < radius)
			point_in_Ball.push_back(root->coordinates);
		return point_in_Ball;
	}
	for(auto child_node : root->children)
		if(utils::vectors_distance(pt,child_node->coordinates)<(radius+child_node->radius))
			for(auto point:pointInBall(child_node, pt, radius))
				point_in_Ball.push_back(point);
	return point_in_Ball;
}