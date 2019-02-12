#include <string>
#include <vector>
#include <set>
#include "simplexArrayList.hpp"
#include "utils.cpp"

// simplexArrayList constructor, currently no needed information for the class constructor
simplexArrayList::simplexArrayList(double maxE, std::vector<std::vector<double>> _distMatrix){
	simplexType = "simplexArrayList";
	maxEpsilon = maxE;
	distMatrix = _distMatrix;
	indexCount = 0;
}

bool sortinrev(const std::pair<double,std::vector<unsigned>> &a, const std::pair<double,std::vector<unsigned>> &b) 
{ 
	return (a.first > b.first); 
} 

double simplexArrayList::getSize(){
	//Calculate size of original data
	
	size_t size = 0;
	
	//Calculate size of edges
	for(int i = 0; i < weightedGraph.size(); i++){
		
		//Size is the ([# weighted graph entries] x [std::pair size]) + ([dimension of graph] * [vector entry size]) 
		size += (weightedGraph[i].size() * sizeof(std::make_pair(0.0,weightedGraph[i][0].second))) + \
					(i * sizeof(unsigned));
	}
	
	//Calculate size of weights
	for(auto row : weights)
		size += sizeof(row);
	std::cout << "SIZE: " << size << std::endl; //<< " : " << sizeof(weightedGraph[0][0].second)*weightedGraph[0][0].second.size() <<std::endl;
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
	
	if(vector.empty())
		return;
	
	//If this is the first point inserted...
	if(weightedGraph.size() == 0){	
		std::vector<unsigned> blank = {0};
		tempGraph = std::make_pair(0.0, blank);
		std::vector<std::pair<double,std::vector<unsigned>>> pair;
		pair.push_back(tempGraph);
		weightedGraph.push_back(pair);
		
		return;
	}
	
	
	//If there are already points, do a brute-force compare
	//		this will take a comparison to every existing point
	else {
		
		auto i = weightedGraph[0].size();
		std::vector<unsigned> blank = {i};
		weightedGraph[0].push_back(std::make_pair(0.0,blank));
		
		//Iterate through each existing to compare to new insertion
		for(unsigned j = 0; j < weightedGraph[0].size(); j++){
			
			auto dist = distMatrix[j][i];
			
			//Filter distances <= maxEpsilon, > 0 (same point)
			if(dist <= maxEpsilon && dist > 0){
				
				//Create an Edge vector (pair) 
				//NOTE: do this in opposite order so pairs are ordered! -> {J, I}
				std::vector<unsigned> edge = {j,i};
				
				if(weightedGraph.size() == 1){
					tempGraph = std::make_pair(dist, edge);
					std::vector<std::pair<double,std::vector<unsigned>>> pair;
					pair.push_back(tempGraph);
					weightedGraph.push_back(pair);
					weights.push_back(dist);
				}
				
				bool testExist = false;
				for(auto z : weightedGraph[1]){
					if(z.second == edge)
						testExist=true;
					}
				if(!testExist)
					weightedGraph[1].push_back(std::make_pair(dist,edge));
				
				if(std::find(weights.begin(), weights.end(), dist) == weights.end())
					weights.push_back(dist);
			}
		}
	}	
	
	return;
}


// Wrapper to expose edges
std::vector<std::vector<unsigned>> simplexArrayList::getEdges(int dim, double epsilon){
	std::vector<std::vector<unsigned>> ret;
	utils ut;
	std::cout << "\t_______" << std::endl;
	std::cout << "\tEpsilon: " << epsilon << "\tdim: " << dim << std::endl;
	
	if(dim - 1 < 0)
		return ret;
	if(weightedGraph.size() < dim)
		return ret;
	
	for(auto a : weightedGraph[dim-1]){
		if(a.first <= epsilon){
			std::cout << "\t";
			ret.push_back(a.second);
			ut.print1DVector(a.second);
		}
	}
	std::cout << "\tRetSize: " << ret.size() << "\tWG Size: " << weightedGraph[dim-1].size() << std::endl;
	
	
	return ret;
}

// Wrapper to expose edges
std::vector<std::pair<double,std::vector<unsigned>>> simplexArrayList::getAllEdges(){
	std::vector<std::pair<double, std::vector<unsigned>>> ret;
	
	for(auto a : weightedGraph){
		for(auto b : a)
			ret.push_back(b);
	}
	
	std::sort(ret.begin(), ret.end(), sortinrev);
	std::cout << "RetSize: " << ret.size() << std::endl;
	return ret;
}

// Search function to find a specific vector in the simplexArrayList
void simplexArrayList::find(std::vector<double> vector){
	
	return;
}

// Output the total simplices stored in the simplical complex
int simplexArrayList::simplexCount(){
	utils ut;
	int simplexRet = 0;
	
	for(auto a : weightedGraph){
		simplexRet += a.size();
	}
	
	return simplexRet;
}

// Output the total vertices stored in the simplical complex
int simplexArrayList::vertexCount(){
	return weightedGraph[0].size();
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
	//Store the complex built from the VR expansion
	std::vector<std::pair<double,std::vector<unsigned>>> weightGraph = weightedGraph[1];
	
	
	//Iterate up to max dimension of simplex, starting at dim 2 (edges)
	for(unsigned i = 2; i <= dim+1; i++){
		if(weightedGraph.size() == i+1)
			weightGraph = weightedGraph[i];
		
		//Store d-dimensional simplices
		std::vector<std::pair<double,std::vector<unsigned>>> test;
		
		//Iterate through each element in the previous dimension's edges
		for(unsigned j = 0; j < weightGraph.size(); j++){
						
			//First search for intersections of the current element
			for(unsigned t = j+1; t < weightGraph.size(); t++){
				auto simp = ut.intersect(weightGraph[j].second, weightGraph[t].second,true);
				

				//This point intersects; potential candidate for a higher-level simplice
				if (simp.first.size() >0){
					//Get the max weight of simplex
					double weight = 0;
					if(weightGraph[j].first > weightGraph[t].first)
						weight = weightGraph[j].first;
					else
						weight = weightGraph[t].first;
					
					for(int k = t+1; k < weightGraph.size(); k++){
						
						if(weightGraph[k].second == simp.first){
							
							if(weightGraph[k].first > weight)
								weight = weightGraph[k].first;
								
							bool testExist = false;
							for(auto z : test){
								if(z.second == simp.second)
									testExist=true;
							}
							if(!testExist)
								test.push_back(std::make_pair(weight,simp.second));
								
							if(std::find(weights.begin(), weights.end(), weight) == weights.end())
								weights.push_back(weight);
						}
					}
				}
			}
		}
		
		//Clear the weightGraph for the next iteration (remove d-1 simplices)
		weightGraph.clear();
		
		//Store the found d-dimensional simplices into weightedGraph and temp (weightGraph)
		for (auto n : test){
			
			if(weightedGraph.size() == i){
				std::vector<std::pair<double,std::vector<unsigned>>> pair;
				pair.push_back(n);
				weightedGraph.push_back(pair);
			}
			else
				weightedGraph[i].push_back(n);
			weightGraph.push_back(n);
		}
		
	} 
	
	
	//Sort the simplices by weight
	for(auto a : weightedGraph){
		std::sort(a.begin(), a.end(), sortinrev);
		std::cout << "************" << std::endl;
		for (auto z : a){
			ut.print1DVector(z.second);
		}
	}
	
	weights.push_back(0.0);
	std::sort(weights.begin(), weights.end(), std::greater<>());
	
	for(auto a : weights){
		std::cout << a << std::endl;
	}
	std::cout << std::endl;
	
	return;
}
