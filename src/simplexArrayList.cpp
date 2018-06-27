#include <string>
#include <vector>
#include "simplexArrayList.hpp"

// simplexArrayList constructor, currently no needed information for the class constructor
simplexArrayList::simplexArrayList(){
	simplexType = "simplexArrayList";
}

double simplexArrayList::getSize(){
	//Calculate size of original data
	
	size_t size = 0;
	
	//Calculate size of edges
	for(auto row : edges){
		for(auto index : row)
			size += sizeof(index);
	}	
	
	//Calculate size of weights
	for(auto row : weights)
		size += sizeof(row);
	
	return size;
}

void simplexArrayList::insert(std::vector<double> vector){
	
	return;
}

void simplexArrayList::find(std::vector<double> vector){
	
	return;
}

int simplexArrayList::simplexCount(){
	
	return -1;
}

int simplexArrayList::vertexCount(){
	
	return -1;
}
