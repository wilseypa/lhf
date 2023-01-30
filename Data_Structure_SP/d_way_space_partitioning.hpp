
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

#include "../Utils/utils.hpp"
#include "../Utils/readInput.hpp"


std::vector<std::vector<double>> referenceHypertetrhedron;
int dim;

class dwaytreenode {
public:
    std::vector<double> coordinates;
    dwaytreenode*parent;
    int parenttochilddirection;
    std::vector<dwaytreenode*> children;
    dwaytreenode()
    {
		this->coordinates.assign(dim, 0);
		this->parent = nullptr;
    }
    dwaytreenode(std::vector<double> coords,int direction)
    {
        this->coordinates = coords;
        this->parenttochilddirection = direction;
    }
    dwaytreenode* buildDwayTree(std::vector<std::vector<double>> data,int direction);
    void printTree(dwaytreenode *root);
    void printLevelOrder(dwaytreenode* root);
    void printCurrentLevel(dwaytreenode* root, int level);
    int height(dwaytreenode* node);
    void pathToCell(std::vector<dwaytreenode*> &path,dwaytreenode* root,dwaytreenode* node);
	std::vector<std::vector<double>>  findCellBoundingPolytope(dwaytreenode* root, dwaytreenode* node);

};

void dwaytreenode:: printLevelOrder(dwaytreenode* root){
    int h = height(root);
    int i;
    for (i = 1; i <= h; i++){
		std::cout<<"Level ::"<<i<<"\n";
        printCurrentLevel(root, i);
	}
}
 
void dwaytreenode:: printCurrentLevel(dwaytreenode* root, int level){
    if (root == nullptr)
        return;
    if (level == 1){
        //for(auto x : root->coordinates)
		//	std::cout<<x<<" ";
		std::cout<<" direction "<<root->parenttochilddirection<<"\n";
	}
    else if (level > 1) {
		for(auto child:root->children){
			printCurrentLevel(child, level - 1);

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

dwaytreenode* dwaytreenode::buildDwayTree(std::vector<std::vector<double>> data,int direction){
	std::vector<double> centroid;
	std::vector<std::vector<double>> coords = utils::transpose(data);
	for(auto coordDim : coords){
		centroid.push_back(utils::getAverage(coordDim));
	}
	dwaytreenode* root = new dwaytreenode(centroid,direction);
	if(data.size()<=1){
		return root;
	}
    //**********************Partition the data into d+1 buckets*******************
    std::vector<double> partitions; 
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
	std::vector<std::vector<std::vector<double>>> dwaypartition(dim+1, std::vector<std::vector<double>>(0, std::vector<double>(0)));
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
			dwaytreenode *child = root->buildDwayTree(part,direction);
			child->parent = root;
			root->children.push_back(child);
		}
		else if(part.size()==1){
			dwaytreenode *child = new dwaytreenode(part[0],direction);
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
		std::cout<<cosine<<" cosine ";
		if(maxvalue<cosine){
			maxvalue = cosine;
			assignedpartition = direction;
		}
		std::cout<<assignedpartition<<" ap";
	}
   std::cout<<" "<<assignedpartition<<" Direction";
   for(auto child : croot->children){
	   if(child->parenttochilddirection==assignedpartition){
		    std::cout<<assignedpartition<<" ap\n";
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
				boundingVertices[x->parenttochilddirection] = x->coordinates;
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
