#include <string>
#include <vector>
#include <set>
#include "simplexArrayList.hpp"
#include "utils.cpp"

// simplexArrayList constructor, currently no needed information for the class constructor
simplexArrayList::simplexArrayList(double maxE){
	simplexType = "simplexArrayList";
	maxEpsilon = maxE;
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
	utils ut;
	
	std::pair<double,std::vector<unsigned>> tempGraph;
	if(data.size() == 0){
		data.push_back(vector);
		return;
	}
	else {
		auto i = data.size() ;
		//Iterate through each existing to compare to new insertion
		for(unsigned j = 0; j < data.size(); j++){
			
				
			//Calculate vector distance 
			auto dist = ut.vectors_distance(vector,data[j]);
			//Filter distances <= epsilon, > 0 (same point)
			if(dist <= maxEpsilon && dist > 0){
				std::vector<unsigned> edge = {i,j};
				edges.push_back(edge);
				weights.push_back(dist);
				tempGraph = std::make_pair(dist,edge);
				weightedGraph.push_back(tempGraph);
			}
		}
		data.push_back(vector);
	}	
	
	//std::cout << "Edges: " << edges.size() << std::endl;
	//std::cout << "Points: " << data.size() << std::endl;
	//std::cout << "Weights: " << weights.size() << std::endl << std::endl;
	
	return;
}

std::vector<std::vector<unsigned>> simplexArrayList::getEdges(int dim, double epsilon){
	return edges;
}

void simplexArrayList::find(std::vector<double> vector){
	
	return;
}

int simplexArrayList::simplexCount(){
	utils ut;
	std::cout << "\tSimplices: " << edges.size() << std::endl;
	
	for(auto q : edges){
		ut.print1DVector(q);
	}
	
	return -1;
}

int simplexArrayList::vertexCount(){
	
	return -1;
}

void simplexArrayList::expandDimensions(int dim){
	utils ut;
	//Store the complex built from the VR expansion
	std::vector<std::vector<unsigned>> loc_edges;
	std::vector<double> loc_weights;
	std::vector<std::vector<unsigned>> loc_temp;
	
	
	for(auto sortTemp : edges){
		sort(sortTemp.begin(), sortTemp.end());
		loc_temp.push_back(sortTemp);
		ut.print1DVector(sortTemp);
	}
	
	std::cout << "\t1-Dimensional Simplices: " << edges.size() << std::endl;
	
	//Iterate up to max dimension of simplex
	for(unsigned i = 2; i < dim; i++){
		std::vector<std::vector<unsigned>> test;
		std::vector<std::pair<std::vector<unsigned>, std::vector<unsigned>>>  tempIntersections;
		
		//Iterate through each element in the previous dimension's edges
		for(unsigned j = 0; j < loc_temp.size(); j++){
			
			//First search for intersections of the current element
			for(unsigned t = j+1; t < loc_temp.size(); t++){
				auto simp = ut.intersect(loc_temp[j], loc_temp[t],false);

				//This point intersects; potential candidate for a higher-level simplice
				if (simp.first.size() >0){//== ix){
					//std::cout<< "pushing simplex" << std::endl;
					tempIntersections.push_back(simp);
				}
			}
			//Loop through the current tempIntersections list
			for(auto tempInt : tempIntersections){				
				
				if(loc_temp[j] == tempInt.first){
					//Copy the tested simplex into the new intersection
					//std::cout<< "Copying Simplex!" << std::endl;
					test.push_back(tempInt.second);
				
				}
			}		
		}		
		//Filter any duplicate simplices
		std::set<std::vector<unsigned>> s(test.begin(), test.end());
		test.assign(s.begin(), s.end());

		for(auto q : test){
			ut.print1DVector(q);
		}
		loc_temp.clear();
		for (auto n : test){
			edges.push_back(n);
			loc_temp.push_back(n);
		}
	}
	
	
	return;
}
