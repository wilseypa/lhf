#include <string>
#include <vector>
#include <set>
#include "simplexArrayList.hpp"

// simplexArrayList constructor, currently no needed information for the class constructor
simplexArrayList::simplexArrayList(double maxE, double maxD, std::vector<std::vector<double>> _distMatrix){
	simplexType = "simplexArrayList";
	maxEpsilon = maxE;
	maxDimension = maxD;
	distMatrix = _distMatrix;
	indexCount = 0;
}

double simplexArrayList::getSize(){
	//Calculate size of original data
	size_t size = 0;
	
	//Calculate size of edges
	for(int i = 0; i < weightedGraph.size(); i++){
		
		//Size is the ([# weighted graph entries] x [std::pair size]) + ([dimension of graph] * [vector entry size]) 
		size += (weightedGraph[i].size() * sizeof(weightedGraph[i][0]));
	}
	
	return size;
}


// Insert for simplexArrayList -> O((n+1)(n+2)/2) -> O(n^2)
//		Sequence: 0 , 1 , 3 , 6 , 10 , 15
//		(AKA very inefficient)
//
void simplexArrayList::insert(std::vector<double> &vector){
	//Create a temporary pair to hold the weight and 2-D edge
	//	e.g.  1.82 , {1, 5} would represent an edge between
	//		points 1 and 5 with a weight of 1.82
		
		
	std::vector<unsigned> vertex = {0};
	if(vector.empty())
		return;
	
	//If this is the first point inserted...
	if(weightedGraph.size() == 0){
		std::vector<std::pair<std::vector<unsigned>,double>> a = {std::make_pair(vertex,0.0)};
		weightedGraph.push_back(a);
		return;
	}
	
	//If there are already points, do a brute-force compare
	//		this will take a comparison to every existing point
	else {
		
		unsigned i = weightedGraph[0].size();
		vertex = {i};
		if(maxDimension > 0){
			//Iterate through each existing to compare to new insertion
			for(unsigned j = 0; j < weightedGraph[0].size(); j++){
				
				auto dist = distMatrix[j][i];
				
				//Filter distances <= maxEpsilon, > 0 (same point)
				if(dist <= maxEpsilon){
					
					//Create an Edge vector (pair) 
					//NOTE: do this in opposite order so pairs are ordered! -> {J, I}
					std::vector<unsigned> edge = {j,i};
					
					if(weightedGraph.size() == 1)
						weightedGraph.push_back({std::make_pair(edge, dist)});
					else if(std::find(weightedGraph[1].begin(), weightedGraph[1].end(), std::make_pair(edge, dist)) == weightedGraph[1].end())
						weightedGraph[1].push_back(std::make_pair(edge,dist));
				}
			}
		}
		weightedGraph[0].push_back(std::make_pair(vertex,0.0));
		
	}	
	
	return;
}


// Wrapper to expose edges
std::vector<std::vector<unsigned>> simplexArrayList::getDimEdges(int dim, double epsilon){
	std::vector<std::vector<unsigned>> ret;
	if(dim - 1 < 0)
		return ret;
	if(weightedGraph.size() < dim)
		return ret;	

	for(auto a : weightedGraph[dim-1]){
		
		bool isTrue = true;
		for(int i = 0; i < a.first.size(); i++){
			for(int j = i+1; j < a.first.size(); j++){
				if(distMatrix[a.first[i]][a.first[j]] > epsilon){
					isTrue =false;
					//weightedGraph[dim-1].erase(weightedGraph[dim-1].begin() + a);
					break;
				}
			}
		}
		
		if(isTrue){
			ret.push_back(a.first);
		}
	}
	
	return ret;
}

// Wrapper to expose edges
std::vector<std::vector<std::pair<std::set<unsigned>, double>>> simplexArrayList::getAllEdges(double epsilon){
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> ret;
	
	for(auto dim : weightedGraph){
		std::vector<std::pair<std::set<unsigned>,double>> dimGraph;
		
		for(auto edge : dim){
			std::set<unsigned> curSet;
			std::copy(edge.first.begin(), edge.first.end(), std::inserter(curSet, curSet.end()));
			dimGraph.push_back(std::make_pair(curSet, edge.second));			
		}
		
		if(dimGraph.size() > 0)
			ret.push_back(dimGraph);
	}
	
	return ret;
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
bool simplexArrayList::find(std::vector<unsigned> vector){
	//Search the weighted graph from the size of the vector
	for(auto v = 0; v < weightedGraph[vector.size() - 1].size(); v++){
		//ut.print1DVector(weightedGraph[vector.size() - 1][v]);
		
		if(weightedGraph[vector.size() - 1][v].first == vector){
			return true;
		}
	}
	return false;
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
double simplexArrayList::findWeight(std::vector<unsigned> vector){
	//Search the weighted graph from the size of the vector
	for(auto v = 0; v < weightedGraph[vector.size() - 1].size(); v++){
		//ut.print1DVector(weightedGraph[vector.size() - 1][v]);
		
		if(weightedGraph[vector.size() - 1][v].first == vector){
			return weightedGraph[vector.size() - 1][v].second;
		}
	}
	return -1;
}

// Output the total simplices stored in the simplical complex
int simplexArrayList::simplexCount(){
	int simplexRet = 0;
	
	for(auto a : weightedGraph){
		simplexRet += a.size();
	}
	
	return simplexRet;
}

// Output the total vertices stored in the simplical complex
int simplexArrayList::vertexCount(){
	if(weightedGraph.size() == 0)
		return 0;
	return weightedGraph[0].size();
}

// Expand the simplexArrayList to incorporate higher-level simplices 
//	-> O(d((n+1)(n+2)/2)) -> O(dn^2) -> where n is the number of d-1 simplices
//		Sequence: 0 , 1 , 3 , 6 , 10 , 15
//		(AKA very inefficient)
//
//	Do this by comparing each simplex to subsequent simplices; if they intersect
//		with a face, search for the remaining faces
//
//
void simplexArrayList::expandDimensions(int dim){	
	
	std::vector<unsigned> tempVect;	
	
	//Iterate up to max dimension of simplex, starting at dim 2 (edges)
	for(unsigned d = 2; d <= dim; d++){
		
		//Check if we need to break from expanding dimensions (no more edges)
		if(weightedGraph.size() < d)
			break;
		
		//Store d-dimensional simplices
		std::vector<std::vector<unsigned>> test;
		
		//Iterate through each element in the current dimension's edges
		for(unsigned j = 0; j < weightedGraph[d-1].size(); j++){
			
			//First search for intersections of the current element
			for(unsigned t = j+1; t < weightedGraph[d-1].size(); t++){
					
				//Symmetric Diff will give us the 
				auto simp = ut.symmetricDiff(weightedGraph[d-1][j].first, weightedGraph[d-1][t].first,true);
				std::vector<unsigned> totalVector = simp;
				double maxWeight = weightedGraph[d-1][j].second > weightedGraph[d-1][t].second ? weightedGraph[d-1][j].second : weightedGraph[d-1][t].second;
				
				
				//This point intersects; potential candidate for a higher-level simplice
				//	Note - this is supposed to be 2, and will always be 2
				if (simp.size() == 2){
					
					bool create = true;
					auto m = ut.setIntersect(weightedGraph[d-1][j].first,weightedGraph[d-1][t].first,true);
					
					//Case that we have a single vertex as the intersect
					if(m.size() == 1){
						std::vector<unsigned> searchVector = simp;
						
						sort(searchVector.begin(), searchVector.end());
						totalVector = searchVector;
						for(auto pt : m)
							totalVector.push_back(pt);
				
						sort(totalVector.begin(), totalVector.end());
						
						double wt = 0;
						if((wt = findWeight(searchVector)) < 0){
							create = false;
						} else if (wt > maxWeight)
							maxWeight = wt;
							
						
					
					} else if(m.size() > 1){
						
						for(auto z : ut.getSubsets(m)){
							
							std::vector<unsigned> searchVector = simp;
							for(unsigned pt : z){
								searchVector.push_back(pt);
							}
							totalVector = ut.setUnion(totalVector, searchVector, true);
							
							double wt = 0;
							
							if((wt = findWeight(searchVector)) < 0){
								create = false;
								break;
							} else if (wt > maxWeight)
								maxWeight = wt;
						}
					}
					
					if(create){
						if(weightedGraph.size() == d){
							std::vector<std::pair<std::vector<unsigned>, double>> tot = {std::make_pair(totalVector, maxWeight)};
							weightedGraph.push_back(tot);
						}else{
							if(std::find(weightedGraph[d].begin(), weightedGraph[d].end(), std::make_pair(totalVector, maxWeight)) == weightedGraph[d].end()){
								weightedGraph[d].push_back(std::make_pair(totalVector, maxWeight));
							}
						}
					}
				}
			}
		}
	} 
	
	//Sort the simplices by weight
	for(auto a : weightedGraph){
		std::sort(a.begin(), a.end(), std::greater<>());
	}
	
	return;
}

void simplexArrayList::reduceComplex(){
	
	//Start with the largest dimension
	ut.writeDebug("simplexArrayList","Reducing complex, starting simplex count: " + std::to_string(simplexCount()));
	
	
	for(auto i = weightedGraph.size()-1; i > 1; i--){
		
		std::vector<std::vector<unsigned>> removals;
		std::vector<std::vector<unsigned>> checked;
		std::vector<unsigned> currentSimplex;

		while(checked.size() != weightedGraph[i].size()){
			for(auto l : weightedGraph[i]){
				if(std::find(checked.begin(),checked.end(),l.first) == checked.end()){
					auto ret = recurseReduce(l, removals, checked);
					removals = ret.first;
					checked = ret.second;
				}
			}
		}	
		
		//Remove the removals
		for(auto rem : removals){
			deletion(rem);
		}
			
	}
	
	ut.writeDebug("simplexArrayList","Finished reducing complex, reduced simplex count: " + std::to_string(simplexCount()));
	
	return;
}

std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<unsigned>>> simplexArrayList::recurseReduce(std::pair<std::vector<unsigned>,double> simplex, std::vector<std::vector<unsigned>> removals, std::vector<std::vector<unsigned>> checked){
	checked.push_back(simplex.first);
	auto subsets = ut.getSubsets(simplex.first);
	std::vector<unsigned> maxFace;
	
	bool canRemove = true;
	
	//Look at each face
	for(auto face : subsets){
		
		//Check if the face is shared; if so, recurse
		for(auto simp : weightedGraph[simplex.first.size() - 1]){
			
			if(simp != simplex && std::find(checked.begin(), checked.end(), simp.first) == checked.end()){
				auto sDiff = ut.symmetricDiff(simp.first, face,true);

				//This point intersects;
				if (sDiff.size() == 1){
					auto ret = recurseReduce(simp, removals, checked);
					removals = ret.first;
					checked = ret.second;
				
					//Check if the simplex was not removed
					if(std::find(removals.begin(), removals.end(), simp.first) == removals.end()){
						canRemove = false;
						break;
					}
				}
				
			}
		}
		
		//Check if the face is the max face
		double wt = -1;
		if((wt = findWeight(face)) == simplex.second){
			maxFace = face;
		}
		
	}	
	
	if(canRemove){
		removals.push_back(simplex.first);
		removals.push_back(maxFace);
	}	
	
	return std::make_pair(removals, checked);
	
}

bool simplexArrayList::deletion(std::vector<unsigned> vector){
	//Search the weighted graph from the size of the vector
	for(auto v = 0; v < weightedGraph[vector.size() - 1].size(); v++){
		//ut.print1DVector(weightedGraph[vector.size() - 1][v]);
		
		if(weightedGraph[vector.size() - 1][v].first == vector){
			weightedGraph[vector.size() - 1].erase(weightedGraph[vector.size() - 1].begin() + v);
			return true;
		}
	}
	return false;
}
