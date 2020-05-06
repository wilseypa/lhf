#include <string>
#include <vector>
#include <set>
#include "simplexArrayList.hpp"

// simplexArrayList constructor, currently no needed information for the class constructor
simplexArrayList::simplexArrayList(double maxE, double maxD, std::vector<std::vector<double>>* _distMatrix){
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
	for(int i = 0; i < simplexList.size(); i++){
		
		//Size is the ([# weighted graph entries] x [std::pair size]) + ([dimension of graph] * [vector entry size]) 
		size += (simplexList[i].size() * sizeof((*simplexList[i].begin())));
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
		
		
	std::set<unsigned> vertex = {0};
	if(vector.empty())
		return;
		
	simplexBase::simplexNode* insNode = new simplexBase::simplexNode();
	
	//If this is the first point inserted...
	if(simplexList.size() == 0){
		insNode->simplex = vertex;
		insNode->weight = 0.0;
		simplexList.push_back({insNode});
		return;
	}
	
	//If there are already points, do a brute-force compare
	//		this will take a comparison to every existing point
	else {
		
		unsigned i = simplexList[0].size();
		vertex = {i};
		if(maxDimension > 0){
			//Iterate through each existing to compare to new insertion
			for(unsigned j = 0; j < simplexList[0].size(); j++){
				
				auto d2 = distMatrix[j];
				auto d3 = d2[i];
				double dist = d3[0];
				
				//Filter distances <= maxEpsilon, > 0 (same point)
				if(dist <= maxEpsilon){
					
					//Create an Edge vector (pair) 
					//NOTE: do this in opposite order so pairs are ordered! -> {J, I}
					std::set<unsigned> edge = {j,i};
					simplexBase::simplexNode* insNode = new simplexBase::simplexNode();
					insNode->simplex = edge;
					insNode->weight = dist;
					
					if(simplexList.size() == 1)
						simplexList.push_back({insNode});
					else
						simplexList[1].insert(insNode);
				}
			}
		}
		
		insNode->simplex = vertex;
		insNode->weight = 0.0;
		simplexList[0].insert(insNode);
		
	}	
	
	return;
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
bool simplexArrayList::find(std::set<unsigned> vector){
	for(auto simplexSetIter = simplexList[vector.size() - 1].begin(); simplexSetIter != simplexList[vector.size() - 1].end(); simplexSetIter++){
		if((*simplexSetIter)->simplex == vector){
			return true;
		}
	}
	return false;
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
double simplexArrayList::findWeight(std::set<unsigned> vector){
	//Search the weighted graph from the size of the vector
	for(auto simplexSetIter = simplexList[vector.size() - 1].begin(); simplexSetIter != simplexList[vector.size() - 1].end(); simplexSetIter++){
		if((*simplexSetIter)->simplex == vector){
			return (*simplexSetIter)->weight;
		}
	}
	return -1;
}

// Output the total simplices stored in the simplical complex
int simplexArrayList::simplexCount(){
	int simplexRet = 0;
	
	for(auto a : simplexList){
		simplexRet += a.size();
	}
	
	return simplexRet;
}

// Output the total vertices stored in the simplical complex
int simplexArrayList::vertexCount(){
	if(simplexList.size() == 0)
		return 0;
	return simplexList[0].size();
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
		if(simplexList.size() < d)
			break;
		
		//Store d-dimensional simplices
		std::vector<std::set<unsigned>> test;
		
		auto jSimplexIter = simplexList[d-1].begin();
		auto tSimplexIter = simplexList[d-1].begin();
		
		//Iterate through each element in the current dimension's edges
		for(unsigned j = 0; j < simplexList[d-1].size(); j++, jSimplexIter++){
			
			//Reset our t iterator to j+1
			tSimplexIter = jSimplexIter;
			tSimplexIter++;
			
			//First search for intersections of the current element
			for(unsigned t = j+1; t < simplexList[d-1].size(); t++, tSimplexIter++){
				
				//Symmetric Diff will give us the 
				auto simp = ut.symmetricDiff((*jSimplexIter)->simplex, (*tSimplexIter)->simplex,true);
				std::set<unsigned> totalVector = simp;
				double maxWeight = (*jSimplexIter)->weight > (*tSimplexIter)->weight ? (*jSimplexIter)->weight : (*tSimplexIter)->weight;
				
				
				//This point intersects; potential candidate for a higher-level simplice
				//	Note - this is supposed to be 2, and will always be 2
				if (simp.size() == 2){
					
					bool create = true;
					auto m = ut.setIntersect((*jSimplexIter)->simplex,(*tSimplexIter)->simplex,true);
					
					//Case that we have a single vertex as the intersect
					if(m.size() == 1){
						std::set<unsigned> searchVector = simp;
						
						totalVector = searchVector;
						for(auto pt : m)
							totalVector.insert(pt);
						
						double wt = 0;
						if((wt = findWeight(searchVector)) < 0){
							create = false;
						} else if (wt > maxWeight)
							maxWeight = wt;
							
						
					
					} else if(m.size() > 1){
						
						for(auto z : ut.getSubsets(m)){
							
							std::set<unsigned> searchVector = simp;
							for(unsigned pt : z){
								searchVector.insert(pt);
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
						simplexNode* tot = new simplexNode();
						tot->simplex = totalVector;
						tot->weight = maxWeight;
						
						if(simplexList.size() == d){
							simplexList.push_back({tot});
						}else{
							simplexList[d].insert(tot);
						}
					}
				}
			}
			
		}
	} 
	
	return;
}

void simplexArrayList::reduceComplex(){
	
	//Start with the largest dimension
	ut.writeDebug("simplexArrayList","Reducing complex, starting simplex count: " + std::to_string(simplexCount()));
	if(simplexList.size() == 0){
		return;
	}
	
	for(auto i = simplexList.size()-1; i > 1; i--){
		
		std::vector<std::set<unsigned>> removals;
		std::vector<std::set<unsigned>> checked;
		std::vector<unsigned> currentSimplex;

		while(checked.size() != simplexList[i].size()){
			for(auto l : simplexList[i]){
				if(std::find(checked.begin(),checked.end(),l->simplex) == checked.end()){
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

std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> simplexArrayList::recurseReduce(simplexNode* simplex, std::vector<std::set<unsigned>> removals, std::vector<std::set<unsigned>> checked){
	checked.push_back(simplex->simplex);
	auto subsets = ut.getSubsets(simplex->simplex);
	std::set<unsigned> maxFace;
	
	bool canRemove = true;
	
	//Look at each face
	for(auto face : subsets){
		
		//Check if the face is shared; if so, recurse
		for(auto simp : simplexList[simplex->simplex.size() - 1]){
			
			if(simp != simplex && std::find(checked.begin(), checked.end(), simp->simplex) == checked.end()){
				auto sDiff = ut.symmetricDiff(simp->simplex, face,true);

				//This point intersects;
				if (sDiff.size() == 1){
					auto ret = recurseReduce(simp, removals, checked);
					removals = ret.first;
					checked = ret.second;
				
					//Check if the simplex was not removed
					if(std::find(removals.begin(), removals.end(), simp->simplex) == removals.end()){
						canRemove = false;
						break;
					}
				}
				
			}
		}
		
		//Check if the face is the max face
		double wt = -1;
		//if((wt = findWeight(face)) == simplex->weight){
		//	maxFace = face;
		//}
		
	}	
	
	if(canRemove){
		removals.push_back(simplex->simplex);
		removals.push_back(maxFace);
	}	
	
	return std::make_pair(removals, checked);
	
}

bool simplexArrayList::deletion(std::set<unsigned> vector){
	//Search the weighted graph from the size of the vector
	for(auto simplexSetIter = simplexList[vector.size() - 1].begin(); simplexSetIter != simplexList[vector.size() - 1].end(); simplexSetIter++){
		if((*simplexSetIter)->simplex == vector){
			simplexList[vector.size() - 1].erase(simplexSetIter);
			return true;
		}
	}
	return false;
}

void simplexArrayList::clear(){
	
	for(auto z : simplexList){
		z.clear();
	}
	
	simplexList.clear();
	
}
