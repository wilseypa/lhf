
/* 
 * Need for data structure for d-way space partioning:
 * Issues with current partioning trees:
 * 1. Existing data structures have worst time complexity in higher dimensions.
 * 2. Data points required explods with dimension 2^d.
 * 3. Have same time complexity as bruit force in higher dimensions.
 * 4. Partioning of space increases with dimensions 2^d
 * 
 * How d-way trees can be fruitful in higher dimensionsl space.
 * 1. The space partioning is leaner to the data diension.
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
    dwaytreenode** children;
    dwaytreenode()
    {
		coordinates.assign(dim, 0);
;
        children = nullptr;
    }
    dwaytreenode(std::vector<double> coords)
    {
        this->coordinates = coords;
        this->children = nullptr;
    }
};
 
class dwaySPTree {
    dwaytreenode* root;
 
public:
    dwaySPTree() { 
		root = nullptr; 
		}
    void initialize(std::vector<std::vector<double>> data);
    void printTree();
};

void dwaySPTree::initialize(std::vector<std::vector<double>> data){
	std::vector<double> centroid;
	std::vector<std::vector<double>> coords = utils::transpose(data);
	for(auto coordDim : coords){
		centroid.push_back(utils::getAverage(coordDim));
	}
	dwaytreenode* newNode = new dwaytreenode(centroid);
    //**********************Partition the data into d+1 buckets*******************
    std::vector<double> partitions; 
	for(auto A : data){
		std::vector<double> A_vec;
		int k =0;
		for(auto a:A){
			A_vec.push_back(a-centroid[k]);
			k++;
		}
		int i = 0;
		int assignedpartion = 0;
		int maxvalue = 0;
		for(auto B : referenceHypertetrhedron){
			double cosine =  utils::cosine_similarity(A_vec,B)*(180/M_PI);
			std::cout<<cosine<<" ";
			if(maxvalue<cosine){
				maxvalue = cosine;
				assignedpartion = i;
			}
			i = i+1;
		}
		partitions.push_back(assignedpartion);
	}
	std::vector<std::vector<std::vector<double>>> dwaypartition(dim+1, std::vector<std::vector<double>>(0, std::vector<double>(0)));
	int i = 0;
	for(auto x : partitions){
		std::cout<<x<<" ";
		dwaypartition[x].push_back(data[i]);
		i = i+1;
	}		

	for(auto x : dwaypartition){
		for(auto y : x){
			for(auto z : y){
				std::cout<<z<<" ";
			}
			std::cout<<std::endl;
		}
		std::cout<<"*******************"<<std::endl;
	}
	for(auto x:centroid){
		std::cout<<x<<" ";
	}
		
    if (root == nullptr) {
        root = newNode;
        return;
    }
    
}
 
void dwaySPTree::printTree(){
	 dwaytreenode* temp = root;
     if (root == nullptr) {
        std::cout << "Tree empty" << std::endl;
        return;
    }
    
}
