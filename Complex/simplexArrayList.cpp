#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "simplexArrayList.hpp"
#include <fstream>

// simplexArrayList constructor, currently no needed information for the class constructor
template<typename nodeType>
simplexArrayList<nodeType>::simplexArrayList(double maxE, double maxD) : bin(0,0) {
	this->simplexType = "simplexArrayList";
	this->maxEpsilon = maxE;
	this->maxDimension = maxD;
}

binomialTable::binomialTable(unsigned n, unsigned k) : v(n+1, std::vector<long long>(k+1, 0)){ //Fast computation of binomials with precomputed table
	v[0][0] = 1;

	for(int i=1; i<=n; i++){
		v[i][0] = 1;
		for(int j=1; j<=k; j++){
			v[i][j] = v[i-1][j-1] + v[i-1][j]; //Pascal's Rule
			if(v[i][j] < 0) throw std::overflow_error("Binomial overflow");
		}
	}
}

long long binomialTable::binom(unsigned n, unsigned k){ //Return binomial coefficient
	if(k>n) return 0;
	return v[n][k];
}

//Hash the set by converting to a remapped index
template<typename nodeType>
long long simplexArrayList<nodeType>::simplexHash(const std::set<unsigned>& simplex){
	long long simplexIndex = 0;
	unsigned i = 0;
	auto it = simplex.begin();
	while(it != simplex.end()){
		simplexIndex += bin.binom(*it - this->simplexOffset, ++i); ///TODO - FIX OFFSET
		if(simplexIndex < 0) throw std::overflow_error("Binomial overflow");
		++it;
	}
	return simplexIndex;
}

template<typename nodeType>
unsigned simplexArrayList<nodeType>::maxVertex(long long simplexHash, unsigned high, unsigned low, unsigned k){
	while(high > low){ //Binary search for the max vertex for this simplex
		unsigned mid = (high + low)/2;
		if(bin.binom(mid, k) <= simplexHash) low = mid + 1;
		else high = mid;
	}
	return high - 1;
}

template<typename nodeType>
std::set<unsigned> simplexArrayList<nodeType>::getVertices(long long simplexHash, int dim, unsigned n){
	std::set<unsigned> v;
	for(unsigned k = dim+1; k>0; k--){ //Get all vertices by repeated binary search for max vertex
		n = maxVertex(simplexHash, n, k-1, k);
		v.insert(n);
		simplexHash -= bin.binom(n, k);
	}
	return v;
}

template<typename nodeType>
void simplexArrayList<nodeType>::prepareCofacets(int dim){
	indexConverter.clear();
	for(auto simplex : this->simplexList[dim+1]){
		indexConverter.insert(std::make_pair(simplex->hash, simplex));
	}
}

template<typename nodeType>
void simplexArrayList<nodeType>::prepareFacets(int dim){
	indexConverter.clear();
	for(auto simplex : this->simplexList[dim-1]){
		indexConverter.insert(std::make_pair(simplex->hash, simplex));
	}
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexArrayList<nodeType>::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<std::shared_ptr<nodeType>, std::shared_ptr<nodeType>>& pivotPairs, bool checkEmergent){
	std::vector<std::shared_ptr<nodeType>> ret;
	int nPts = this->simplexList[0].size();
	unsigned k = simplex.size() + 1;
	std::set<unsigned>::reverse_iterator it = simplex.rbegin();
	long long index = simplexHash(simplex);

	//Try inserting other vertices into the simplex
	for(unsigned i=nPts; i-- != 0; ){
		if(it != simplex.rend() && i == *it - this->simplexOffset){ //Vertex i is already in the simplex
			//Now adding vertices less than i -> i is now the kth largest vertex in the simplex instead of the (k-1)th
			index -= bin.binom(i, k-1);
			index += bin.binom(i, k); //Recompute the index accordingly

			--k;	//Now need to check for the (k-1)th vertex in the simplex
			++it; 	//Check for the previous vertex in the simplex (it is a reverse iterator)
		} else{
			auto tempNode = indexConverter.find(index + bin.binom(i, k));
			if(tempNode != indexConverter.end()){ //If this is a valid simplex, add it to the heap
				ret.push_back(tempNode->second);

				//If we haven't found an emergent candidate and the weight of the maximal cofacet is equal to the simplex's weight
				//		we have identified an emergent pair; at this point we can break because the interval is born and dies at the
				//		same epsilon
				if(checkEmergent && tempNode->second->weight == simplexWeight){
					//Check to make sure the identified cofacet isn't a pivot
					if(pivotPairs.find(tempNode->second) == pivotPairs.end()) return ret;
					checkEmergent = false;
				}
			}
		}
	}

	return ret;
}

template<typename nodeType> 
std::vector<nodeType*> simplexArrayList<nodeType>::getAllCofacets(std::shared_ptr<nodeType> simp, const std::unordered_map<long long, std::shared_ptr<nodeType>>& pivotPairs, bool checkEmergent, bool recordVertices, unsigned dim){
	//Method builds out cofacets for incrementalPersistence

	std::vector<nodeType*> ret;
	std::set<unsigned> vertices;

	if(recordVertices) vertices = simp->simplex;
	else vertices = getVertices(simp->hash, dim, this->simplexList[0].size());

	unsigned k = vertices.size() + 1;
	std::set<unsigned>::reverse_iterator it = vertices.rbegin();
	long long index = simp->hash;

	// Try inserting other vertices into the simplex
	for(unsigned i = this->simplexList[0].size(); i-- != 0; ){
		if(it != vertices.rend() && i == *it - this->simplexOffset){ //Vertex i is already in the simplex
			//Now adding vertices less than i -> i is now the kth largest vertex in the simplex instead of the (k-1)th
			index -= bin.binom(i, k-1);
			index += bin.binom(i, k); //Recompute the index accordingly
			--k;	//Now need to check for the (k-1)th vertex in the simplex
			++it; 	//Check for the previous vertex in the simplex (it is a reverse iterator)
		} else{
			double maxWeight = simp->weight;
			for(auto pt : vertices){
				maxWeight = std::max(maxWeight, (*this->distMatrix)[std::min(i, pt)][std::max(i, pt)]);
			}

			if(maxWeight <= this->maxEpsilon){ //Valid simplex
				nodeType* x = new nodeType();
				if(recordVertices){
					x->simplex = vertices;
					x->simplex.insert(i);
				}
				x->weight = maxWeight;
				x->hash = index + bin.binom(i, k);
				ret.push_back(x);

				if(checkEmergent && maxWeight == simp->weight){
					if(pivotPairs.find(index + bin.binom(i, k)) == pivotPairs.end()) return ret;
					checkEmergent = false;
				}
			}
		}
	}

	return ret;
}

template<typename nodeType>
std::vector<nodeType*> simplexArrayList<nodeType>::getAllCofacets(std::shared_ptr<nodeType> simp){
	return getAllCofacets(simp, std::unordered_map<long long, std::shared_ptr<nodeType>>(), false, true, 0);
}

template<typename nodeType>
std::vector<nodeType*> simplexArrayList<nodeType>::getAllFacets(nodeType* simp, bool recordVertices, unsigned dim){
	std::vector<nodeType*> ret;
	std::set<unsigned> vertices;

	if(recordVertices) vertices = simp->simplex;
	else vertices = getVertices(simp->hash, dim + 1, this->simplexList[0].size());

	long long index = simp->hash;
	unsigned k = vertices.size();

	for(auto it = vertices.rbegin(); it != vertices.rend(); ++it){
		unsigned pt = *it;
		double maxWeight = 0;

		for(auto it = vertices.begin(); it != vertices.end(); ++it){
			if(*it == pt) continue;
			for(auto it2 = it; ++it2 != vertices.end(); ){
				if(*it2 == pt) continue;
				maxWeight = std::max(maxWeight, (*this->distMatrix)[*it][*it2]);
			}
		}

		nodeType* x = new nodeType();
		x->weight = maxWeight;

		if(recordVertices){
			x->simplex = vertices;
			x->simplex.erase(x->simplex.find(pt));
		}

		index -= bin.binom(pt, k);
		x->hash = index;
		index += bin.binom(pt, --k);

		ret.push_back(x);
	}

	return ret;
}

template<typename nodeType>
std::vector<nodeType*> simplexArrayList<nodeType>::getAllFacets(std::shared_ptr<nodeType> simp, bool recordVertices, unsigned dimension){
	return getAllFacets(simp.get(), recordVertices, dimension);
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexArrayList<nodeType>::getAllFacets_P(std::shared_ptr<nodeType> simp){
	std::vector<std::shared_ptr<nodeType>> ret;

	long long index = simp->hash;
	unsigned k = simp->simplex.size();

	for(auto it = simp->simplex.rbegin(); it != simp->simplex.rend(); ++it){
		unsigned pt = *it;

		index -= bin.binom(pt, k);

		auto tempNode = indexConverter.find(index);
		if(tempNode != indexConverter.end()){ //If this is a valid simplex, add it to the heap
			ret.push_back(tempNode->second);
		}

		index += bin.binom(pt, --k);
	}

	return ret;
}

template<typename nodeType>
double simplexArrayList<nodeType>::getSize(){
	//Calculate size of original data
	size_t size = 0;

	//Calculate size of edges
	for(int i = 0; i < this->simplexList.size(); i++){

		//Size is the ([# weighted graph entries] x [std::pair size]) + ([dimension of graph] * [vector entry size])
		size += (this->simplexList[i].size() * sizeof((*this->simplexList[i].begin())));
	}

	return size;
}


// Insert for simplexArrayList -> O((n+1)(n+2)/2) -> O(n^2)
//		Sequence: 0 , 1 , 3 , 6 , 10 , 15
//
template<typename nodeType>
void simplexArrayList<nodeType>::insert(){
	//If this is the first point inserted...
	if(this->simplexList.size() == 0) this->simplexList.push_back({});

	unsigned i = this->simplexList[0].size();

	std::shared_ptr<nodeType> insNode = std::make_shared<nodeType>(nodeType({i}, 0.0));
	insNode->hash = i;
	this->simplexList[0].insert(insNode);
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
template<typename nodeType>
bool simplexArrayList<nodeType>::find(std::set<unsigned> vector){
	if (this->simplexList.size() == 0) return false;

	for(auto simplexSetIter = this->simplexList[vector.size() - 1].begin(); simplexSetIter != this->simplexList[vector.size() - 1].end(); simplexSetIter++){
		if((*simplexSetIter)->simplex == vector){
			return true;
		}
	}
	return false;
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
template<typename nodeType>
double simplexArrayList<nodeType>::findWeight(std::set<unsigned> vector){
	//Search the weighted graph from the size of the vector
	for(auto simplexSetIter = this->simplexList[vector.size() - 1].begin(); simplexSetIter != this->simplexList[vector.size() - 1].end(); simplexSetIter++){
		if((*simplexSetIter)->simplex == vector){
			return (*simplexSetIter)->weight;
		}
	}
	return -1;
}

// Output the total simplices stored in the simplical complex
template<typename nodeType>
int simplexArrayList<nodeType>::simplexCount(){
	int simplexRet = 0;

	for(auto a : this->simplexList){
		simplexRet += a.size();
	}

	return simplexRet;
}

// Output the total vertices stored in the simplical complex
template<typename nodeType>
int simplexArrayList<nodeType>::vertexCount(){
	if(this->simplexList.size() == 0)
		return 0;
	return this->simplexList[0].size();
}

template<typename nodeType>
void simplexArrayList<nodeType>::initBinom(){
	bin = binomialTable(this->simplexList[0].size(), this->maxDimension+1);
}

// Expand the simplexArrayList to incorporate higher-level simplices
//	-> O(dnk) -> where n is the number of points, d is the dimension, and k is the number of d-1 simplices
//
//	Do this by comparing each simplex with points to insert
//
template<typename nodeType>
void simplexArrayList<nodeType>::expandDimensions(int dim){
	initBinom();

	//Iterate up to max dimension of simplex, starting at dim 2 (edges)
	for(unsigned d = 1; d <= dim; d++){

		//Check if we need to break from expanding dimensions (no more edges)
		if(this->simplexList.size() < d) break;
		if(this->simplexList.size() == d) this->simplexList.push_back({});

		//Iterate through each element in the current dimension's edges
		for(auto it = this->simplexList[d-1].begin(); it != this->simplexList[d-1].end(); it++){
			//Iterate over points to possibly add to the simplex
			//Use points larger than the maximal vertex in the simplex to prevent double counting
			unsigned minPt = *(*it)->simplex.rbegin() + 1;

			for(unsigned pt = minPt; pt < this->simplexList[0].size(); pt++){
				//Compute the weight using all edges
				double maxWeight = (*it)->weight;
				for(auto i : (*it)->simplex) maxWeight = std::max(maxWeight, (*this->distMatrix)[i][pt]);
			
//***************************For beta complex valid simplex Condition ****************************			
            bool valid = true;
            if(this->complexType == "alphaComplex"){
                  for(auto i : (*it)->simplex) 
						if(!(*this->incidenceMatrix)[i][pt]){
						   valid = false;
						   break;
					   }
            if(valid && maxWeight <= this->maxEpsilon){
					std::shared_ptr<nodeType> tot = std::make_shared<nodeType>(nodeType((*it)->simplex, maxWeight));
					tot->simplex.insert(pt);
					tot->hash = (*it)->hash + bin.binom(pt, tot->simplex.size());
					this->simplexList[d].insert(tot);
			}
//************************************************************************************************
			}else if(maxWeight <= this->maxEpsilon){ //Valid simplex
					std::shared_ptr<nodeType> tot = std::make_shared<nodeType>(nodeType((*it)->simplex, maxWeight));
					tot->simplex.insert(pt);
					tot->hash = (*it)->hash + bin.binom(pt, tot->simplex.size());
					this->simplexList[d].insert(tot);
				}
			}
		}
	}
std::ofstream out("incedenceMatrix2DBeta0.9.csv");
int di=0;
for( auto x : this->simplexList)
	std::cout<<"Count of "<<di++<<"-simplex ::"<<x.size()<<"\n";

std::vector<std::vector<unsigned>> edges(this->simplexList[0].size(),std::vector<unsigned>(this->simplexList[0].size(),0));
for(auto x : this->simplexList[1]){
	std::vector<unsigned> edge;
	for(auto y : x->simplex)
		edge.push_back(y);
    edges[edge[0]][edge[1]] = 1;
	}


for (auto row : edges) {
  for (auto col : row)
    out << col <<' ';
  out << '\n';
}
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexArrayList<nodeType>::expandDimension(std::vector<std::shared_ptr<nodeType>> edges, bool recordVertices, unsigned dim){
	std::vector<std::shared_ptr<nodeType>> nextEdges;
	//Iterate through each element in the current dimension's edges
	for(auto it = edges.begin(); it != edges.end(); it++){
		std::set<unsigned> vertices;

		if(recordVertices) vertices = (*it)->simplex;
		else vertices = getVertices((*it)->hash, dim - 1, this->simplexList[0].size());

		//Iterate over points to possibly add to the simplex
		//Use points larger than the maximal vertex in the simplex to prevent double counting
		unsigned minPt = *vertices.rbegin() + 1;

		for(unsigned pt = minPt; pt < this->simplexList[0].size(); pt++){
			//Compute the weight using all edges
			double maxWeight = (*it)->weight;
			for(auto i : vertices) maxWeight = std::max(maxWeight, (*this->distMatrix)[i][pt]);
			
//***************************For beta complex valid simplex Condition ****************************			
            bool valid = true;
            if(this->complexType == "alphaComplex"){
                  for(auto i : vertices) 
						if(!(*this->incidenceMatrix)[i][pt]){
						   valid = false;
						   break;
					   }
            if(valid){
				std::shared_ptr<nodeType> tot = std::make_shared<nodeType>(nodeType());
				if(recordVertices){
					tot->simplex = vertices;
					tot->simplex.insert(pt);
				}
				tot->weight = maxWeight;
				tot->hash = (*it)->hash + bin.binom(pt, (recordVertices ? tot->simplex.size() : dim + 1));
				nextEdges.push_back(tot);
			}
//************************************************************************************************
			}else if(maxWeight <= this->maxEpsilon){ //Valid simplex
				std::shared_ptr<nodeType> tot = std::make_shared<nodeType>(nodeType());
				if(recordVertices){
					tot->simplex = vertices;
					tot->simplex.insert(pt);
				}
				tot->weight = maxWeight;
				tot->hash = (*it)->hash + bin.binom(pt, (recordVertices ? tot->simplex.size() : dim + 1));
				nextEdges.push_back(tot);
			}
		}
	}

	if(recordVertices) std::sort(nextEdges.begin(), nextEdges.end(), cmpByWeight<std::shared_ptr<nodeType>>());
	return nextEdges;
}

template<typename nodeType>
bool simplexArrayList<nodeType>::deletion(std::set<unsigned> vector){
	//Search the weighted graph from the size of the vector
	for(auto simplexSetIter = this->simplexList[vector.size() - 1].begin(); simplexSetIter != this->simplexList[vector.size() - 1].end(); simplexSetIter++){
		if((*simplexSetIter)->simplex == vector){
			this->simplexList[vector.size() - 1].erase(simplexSetIter);
			return true;
		}
	}
	return false;
}

template<typename nodeType>
simplexArrayList<nodeType>::~simplexArrayList(){
	this->simplexList.clear();
}

//Explicit Template Class Instantiation
template class simplexArrayList<simplexNode>;
template class simplexArrayList<alphaNode>;
template class simplexArrayList<witnessNode>;
