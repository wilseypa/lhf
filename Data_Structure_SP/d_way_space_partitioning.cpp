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
	std::cout<<"Enter Input File Name :: ";
	std::string filename;
	std::getline (std::cin,filename);
	data = rs.readCSV(filename);
    dim =data[0].size();
	referenceHypertetrhedron = utils::genCoordsRegularSimplex(dim);
 
    dwaytreenode *tree; 
    tree = tree->buildDwayTree(data,-1);
    //tree->printTree(tree);
	tree->printLevelOrder(tree);

    return 0;
}
