#include <string>
#include <vector>
#include "pipePacket.hpp"

// pipePacket constructor, currently no needed information for the class constructor
pipePacket::pipePacket(){

}

double pipePacket::getSize(){
	
	//Calculate size of original data
	
	size_t size = 0;
	
	for(auto row : workData.originalData){
		for(auto index : row)
			size += sizeof(index);
	}
	
	//Calculate size of working data
	for(auto row: workData.workingData){
		for(auto index : row)
			size += sizeof(index);
	}
	
	//Calculate size of edges
	for(auto row : workData.edges){
		for(auto index : row)
			size += sizeof(index);
	}	
	
	//Calculate size of weights
	for(auto row : workData.weights)
		size += sizeof(row);
	
	
	return size;
}
