#include <string>
#include <vector>
#include <set>
#include "simplexArrayList.hpp"

// simplexArrayList constructor, currently no needed information for the class constructor
simplexArrayList::simplexArrayList(double maxE, std::vector<std::vector<double>> _distMatrix){
	simplexType = "simplexArrayList";
	maxEpsilon = maxE;
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
	//
	//
	
	if(vector.empty())
		return;
	
	//If this is the first point inserted...
	if(weightedGraph.size() == 0){	
		weightedGraph.push_back({{0}});
		return;
	}
	
	//If there are already points, do a brute-force compare
	//		this will take a comparison to every existing point
	else {
		
		auto i = weightedGraph[0].size();
		
		//Iterate through each existing to compare to new insertion
		for(unsigned j = 0; j < weightedGraph[0].size(); j++){
			
			auto dist = distMatrix[j][i];
			
			//Filter distances <= maxEpsilon, > 0 (same point)
			if(dist <= maxEpsilon){
				
				//Create an Edge vector (pair) 
				//NOTE: do this in opposite order so pairs are ordered! -> {J, I}
				std::vector<unsigned> edge = {j,i};
				
				if(weightedGraph.size() == 1)
					weightedGraph.push_back({edge});
				else if(std::find(weightedGraph[1].begin(), weightedGraph[1].end(), edge) == weightedGraph[1].end())
					weightedGraph[1].push_back(edge);
			}
		}
		weightedGraph[0].push_back({i});
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
		for(int i = 0; i < a.size(); i++){
			for(int j = i+1; j < a.size(); j++){
				if(distMatrix[a[i]][a[j]] > epsilon){
					isTrue =false;
					//weightedGraph[dim-1].erase(weightedGraph[dim-1].begin() + a);
					break;
				}
			}
		}
		
		if(isTrue){
			ret.push_back(a);
		}
	}
	
	return ret;
}

// Wrapper to expose edges
std::vector<std::vector<std::pair<std::set<unsigned>, double>>> simplexArrayList::getAllEdges(double epsilon){
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> ret;
	
	/*for(int dim = 0; dim < weightedGraph.size(); dim++){
		
		std::vector<std::pair<std::vector<unsigned>, double>> temp;
		for(int a = 0; a < weightedGraph[dim].size(); a++){
			if(weightedGraph[dim][a].size() > 0){ 
				double maxDist = 0;
				bool isTrue = true;
				for(int t = 0; t < weightedGraph[dim][a].size(); t++){
					for(int s = t+1; s < weightedGraph[dim][a].size(); s++){
						
						double curDist = distMatrix[weightedGraph[dim][a][t]][weightedGraph[dim][a][s]];
						
						if(curDist > epsilon){
							
							isTrue =false;
							//weightedGraph[dim].erase(weightedGraph[dim].begin() + a);
							//weightedGraph[dim].erase(remove(weightedGraph[dim].begin(), weightedGraph[dim].end(), weightedGraph[dim][a]));
							break;
						} else if (curDist > maxDist)
							maxDist = curDist;
					}
					if(!isTrue){
						break;
					}
				}
				
				if(isTrue){
					temp.push_back(std::make_pair(weightedGraph[dim][a], maxDist));
				}
			}
		}
		ret.push_back(temp);
	}*/
	
	
	return ret;
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
bool simplexArrayList::find(std::vector<unsigned> vector){
	
	for(auto d = 0; d = weightedGraph.size()-1; d++){
		for(auto v = 0; v = weightedGraph[d].size(); v++){
			for(auto p = 0; p = weightedGraph[d][v].size(); p++){
				if(weightedGraph[d][v][p] = vector[p]){
					return true;
				}

			}

		}

	}
	return false;
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
//
void simplexArrayList::expandDimensions(int dim){
	std::vector<unsigned> tempVect;	
	
	
	//Iterate up to max dimension of simplex, starting at dim 2 (edges)
	for(unsigned d = 2; d <= dim; d++){
		
		//Store d-dimensional simplices
		std::vector<std::vector<unsigned>> test;
		
		//Iterate through each element in the previous dimension's edges
		for(unsigned j = 0; j < weightedGraph[d-1].size(); j++){
			//First search for intersections of the current element
			for(unsigned t = j+1; t < weightedGraph[d-1].size(); t++){
				
				auto simp = ut.symmetricDiff(weightedGraph[d-1][j], weightedGraph[d-1][t],true);
				//This point intersects; potential candidate for a higher-level simplice
				//std::cout << "d";
				if (simp.size() == d){
					//std::cout << "t" << std::endl;
					for(int k = t+1; k < weightedGraph[d-1].size(); k++){
						
						if(weightedGraph[d-1][k] == simp){
							simp = ut.setUnion(weightedGraph[d-1][j], weightedGraph[d-1][t],true);
							
							if(std::find(weightedGraph[d].begin(), weightedGraph[d].end(), simp) == weightedGraph[d].end()){
								if(weightedGraph.size() == d)
									weightedGraph.push_back({simp});
								else
									weightedGraph[d].push_back(simp);
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
	
	return;
}


