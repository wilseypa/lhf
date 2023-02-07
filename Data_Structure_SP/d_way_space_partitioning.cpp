/* 
 * Contains Suporting Functions for d-way space partioning Data Structure.
 * 
 */
#include <iostream>
#include "d_way_space_partitioning.hpp"
#include <string>
int main(){
	auto rs = readInput();
	std::vector<std::vector<double>> data;
	//std::cout<<"Enter Input File Name :: ";
	std::string filename = "inputfile3.txt";
	//std::getline (std::cin,filename);
	data = rs.readCSV(filename);
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
    tree = tree->buildDwayTree(data,-1);
    //tree->printTree(tree);
	tree->printLevelOrder(tree,tree);
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
