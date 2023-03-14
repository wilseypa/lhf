/* 
 * Contains Suporting Functions for d-way space partioning Data Structure.
 * 
 */
#include <iostream>
#include "d_way_space_partitioning.hpp"
#include <string>
#include <bits/stdc++.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
using namespace std::chrono;
void fill_row(std::vector<double> & row)
{
    std::generate(row.begin(), row.end(), [](){ return (rand() % 1000000)/1000000; }); 
}

void fill_matrix(std::vector<std::vector<double>> & mat)
{
    std::for_each(mat.begin(), mat.end(), fill_row);
}
    
int main(){
	auto rs = readInput();
	std::vector<std::vector<double>> data;
	//std::cout<<"Enter Input File Name :: ";
	//std::string filename = "filenew.txt";
	std::string filename = "inputfile.txt";
	//std::getline (std::cin,filename);
	data = rs.readCSV(filename);
	int count, dimension;
    //std::cout<<"Enter Size and Dimension of Original Data"<<std::endl;
    //std::cout<<"Size ::";
    //std::cin>>count;
    //std::cout<<"dimension ::";
    //std::cin>>dimension;
    //std::vector<std::vector<double>> data(count, std::vector<double>(dimension, 0));
    //fill_matrix(data);
    /*for (int i=0; i<count; i++){
        for (int j=0; j<dimension; j++)
			std::cout<<data[i][j]<<" ";
		std::cout<<std::endl;
    }
    */
    /*std::cout<<"Enter Size of test Data"<<std::endl;
    std::cout<<"Size ::";
    std::cin>>count;
    std::vector<std::vector<double>> testdata(count, std::vector<double>(dimension, 0));
    fill_matrix(testdata);
    */
    /*for (int i=0; i<count; i++){
        for (int j=0; j<dimension; j++)
			std::cout<<testdata[i][j]<<" ";
		std::cout<<std::endl;
    }
    * */
    dim =data[0].size();	
	referenceHypertetrhedron = utils::genCoordsRegularSimplex(dim,100);
    /*for(auto x:referenceHypertetrhedron){
		for(auto y:x){
			std::cout<<y<<" ";
		}
		std::cout<<"\n";
	}
	*/
	dwaytreenode *tree; 
    tree = tree->buildDwayTree(data,-1,"nary",2);
    //tree->printTree(tree);
   	//tree->printLevelOrder(tree,tree);
   	int hdim;
   	std::cout<<"Enter Homology Dimension";
   	std::cin>>hdim;
   	auto mesh = tree->meshGeneration(tree,tree,1,hdim);
/*    std::cout<<"simplices::\n";
    for(auto x:mesh.first){
		for(auto y:x){
			std::cout<<" ";
			for(auto z:y){
				std::cout<<" "<<z;
			}
			std::cout<<std::endl;
		}
		std::cout<<std::endl;
	}
	std::cout<<"Points\n";
	for(auto x:mesh.second){
		for(auto y:x){
			std::cout<<y<<" ";
		}
		std::cout<<std::endl;
	}
*/
	/*
	std::cout<<tree->checkPointInBall(tree,{0,0} ,1.1,{data[0]});
	auto temp=tree->pointInBall(tree,{0,0},1.1);
	for(auto i:temp){
		for(auto j:i)
			std::cout<<j<<" ";
		std::cout<<std::endl;
	} */

	/*std::vector<double> pt;
	double k;
	std::cin>>k;
	pt.push_back(k);
	std::cin>>k;
	pt.push_back(k);
	*/
	/*
	kdTree kdtree(data, data.size()); //KDTree for efficient nearest neighbor search
    double timekdtree1=0;
    double timedwaytree1=0;
    
	for(auto pt:testdata){
		auto start = high_resolution_clock::now();
		dwaytreenode* NN = tree->findNearestNeighbor(tree, pt);
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
		timedwaytree1 +=duration.count();
		//std::cout<<NN->coordinates[0]<<" "<<NN->coordinates[1]<<" "<<NN->radius;
		auto start1 = high_resolution_clock::now();
		std::vector<double> x = kdtree.nearestPoint(pt);
		auto stop1 = high_resolution_clock::now();
		auto duration1 = duration_cast<microseconds>(stop1 - start1);
		timekdtree1 +=duration1.count();
		
		//std::cout<<x[0]<<" "<<x[1]<<" ";
		if(NN->coordinates==x)
		    std::cout<<"True";
		else
		    std::cout<<"False";
	}
	std::cout<<"\nTime taken by kdtree is  "<<timekdtree1;
	
	std::cout<<"\nTime taken by dway tree is  "<<timedwaytree1;
```*/
	/*std::vector<dwaytreenode*> path;
	dwaytreenode* node = tree->children[0]->children[2]->children[2]->children[0]->children[1];
    //tree->pathToCell(path,tree, node);
    auto boundingVertices =  tree-> findCellBoundingPolytope(tree,node);
    for(auto x:boundingVertices){
		  for(auto z:x){
		     std::cout<<z<<" ";
		  }
		  std::cout<<"\n";
	 }
	 
	std::vector<std::vector<double>> pts {{1, 2, -2},{3, -2, 1},{5, 1, -4}};
	std::vector<double> interior{10,10,10};
	auto hyp = utils::generateHyperplaneFromVertices( pts,interior);
    for(auto x:hyp.first)
		std::cout<<x<<" ";
	std::cout<<hyp.second;
	std::cout<<"Eigen\n";
    int k=0;
    std::cin>>k;
    for(auto x: path){
		std::cout<<x->parenttochilddirection<<" ";
	}*/
    return 0;
}
