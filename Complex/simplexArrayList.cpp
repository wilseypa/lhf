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



// Insert for simplexArrayList -> O((n+1)(n+2)/2) -> O(n^2)
//		Sequence: 0 , 1 , 3 , 6 , 10 , 15
//		(AKA very inefficient)
//
void simplexArrayList::insert(std::vector<double> vector){
	utils ut;
	
	//Create a temporary pair to hold the weight and 2-D edge
	//	e.g.  1.82 , {1, 5} would represent an edge between
	//		points 1 and 5 with a weight of 1.82
	//
	//
	std::pair<double,std::vector<unsigned>> tempGraph;
	
	//If this is the first point inserted...
	if(data.size() == 0){
		data.push_back(vector);
		return;
	}
	//If there are already points, do a brute-force compare
	//		this will take a comparison to every existing point
	else {
		
		auto i = data.size() ;
		
		//Iterate through each existing to compare to new insertion
		for(unsigned j = 0; j < data.size(); j++){
			
			//Calculate vector distance 
			auto dist = ut.vectors_distance(vector,data[j]);
			
			//Filter distances <= maxEpsilon, > 0 (same point)
			if(dist <= maxEpsilon && dist > 0){
				
				//Create an Edge vector (pair) 
				//NOTE: do this in opposite order so pairs are ordered! -> {J, I}
				std::vector<unsigned> edge = {j,i};
				
				//Push edges, weights (deprecated)
				edges.push_back(edge);
				weights.push_back(dist);
				
				//Push edges, weights to weighted graph
				weightedGraph.push_back(std::make_pair(dist,edge));
			}
		}
		
		//After we're finished comparing to existing vectors,
		//	push back the new vector into our stored data
		data.push_back(vector);
	}	
	
	return;
}


// Wrapper to expose edges
std::vector<std::pair<double,std::vector<unsigned>>> simplexArrayList::getEdges(int dim, double epsilon){
	return weightedGraph;
}

// Search function to find a specific vector in the simplexArrayList
void simplexArrayList::find(std::vector<double> vector){
	
	return;
}

// Output the total simplices stored in the simplical complex
int simplexArrayList::simplexCount(){
	utils ut;
	
	return weightedGraph.size();
}

// Output the total vertices stored in the simplical complex
int simplexArrayList::vertexCount(){
	return data.size();
}

// Expand the simplexArrayList to incorporate higher-level simplices 
//	-> O(d((n+1)(n+2)/2)) -> O(dn^2) -> where n is the number of d-1 simplices
//		Sequence: 0 , 1 , 3 , 6 , 10 , 15
//		(AKA very inefficient)
//
//
void simplexArrayList::expandDimensions(int dim){
	utils ut;

	std::vector<unsigned> tempVect;
	//Add the original vectors to the weighted graph
	for(auto i = 0; i < data.size(); i++){
		weightedGraph.push_back(std::make_pair(0.0,tempVect={i}));
	}

	//Store the complex built from the VR expansion
	std::vector<std::pair<double,std::vector<unsigned>>> weightGraph = weightedGraph;
	
	//Sort the 2D simplices by weight
	std::sort(weightedGraph.begin(), weightedGraph.end());
	
	std::cout << "\t1-Dimensional Simplices: " << weightGraph.size() << std::endl;
	
	//Iterate up to max dimension of simplex, starting at dim 2 (edges)
	for(unsigned i = 2; i < dim; i++){
		
		//Store d-dimensional simplices
		std::vector<std::pair<double,std::vector<unsigned>>> test;
		
		//Store vector intersections to compare to as we go
		std::vector<std::pair<double,std::pair<std::vector<unsigned>, std::vector<unsigned>>>>  tempIntersections;
		
		//Iterate through each element in the previous dimension's edges
		for(unsigned j = 0; j < weightGraph.size(); j++){
			
			//First search for intersections of the current element
			for(unsigned t = j+1; t < weightGraph.size(); t++){
				auto simp = ut.intersect(weightGraph[j].second, weightGraph[t].second,false);

				//This point intersects; potential candidate for a higher-level simplice
				if (simp.first.size() >0){
					//Get the max weight of simplex
					double weight = 0;
					
					if(weightGraph[j].first > weightGraph[t].first)
						weight = weightGraph[j].first;
					else
						weight = weightGraph[t].first;
					
					auto weightPair = std::make_pair(weight,simp);
					if(std::find(tempIntersections.begin(), tempIntersections.end(), weightPair) != tempIntersections.end())		
						tempIntersections.push_back(weightPair);
				}
			}
			//Loop through the current tempIntersections list
			for(auto tempInt : tempIntersections){				
				
				if(weightGraph[j].second == tempInt.second.first){
					double weight = 0;
					
					//Copy the tested simplex into the new intersection
					if(weightGraph[j].first > tempInt.first)
						weight = weightGraph[j].first;
					else
						weight = tempInt.first;
					
					test.push_back(std::make_pair(weight,tempInt.second.second));
				}
			}		
		}
		
		//Clear the weightGraph for the next iteration (remove d-1 simplices)
		weightGraph.clear();
		
		//Sort the entries by weight
		std::sort(test.begin(), test.end());
		
		std::cout << "WEIGHT GRAPH" << std::endl << std::endl;
		//Print out the vectors (for clarity)
		for(auto q : test){
			std::cout << q.first << "\t";
			ut.print1DVector(q.second);
		}
		std::cout << std::endl;
		
		//Store the found d-dimensional simplices into weightedGraph and temp (weightGraph)
		for (auto n : test){
			weightedGraph.push_back(n);
			weightGraph.push_back(n);
		}
	}
	
	return;
}
