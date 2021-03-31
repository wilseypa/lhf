#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "simplexArrayList.hpp"

// simplexArrayList constructor, currently no needed information for the class constructor
simplexArrayList::simplexArrayList(double maxE, double maxD) : bin(0,0) {
	simplexType = "simplexArrayList";
	maxEpsilon = maxE;
	maxDimension = maxD;
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
long long simplexArrayList::simplexHash(const std::set<unsigned>& simplex){
	long long simplexIndex = 0;
	unsigned i = 0;
	auto it = simplex.begin();
	while(it != simplex.end()){
		simplexIndex += bin.binom(*it - simplexOffset, ++i); ///TODO - FIX OFFSET
		if(simplexIndex < 0) throw std::overflow_error("Binomial overflow");
		++it;
	}
	return simplexIndex;
}

unsigned simplexArrayList::maxVertex(long long simplexHash, unsigned high, unsigned low, unsigned k){
	while(high > low){ //Binary search for the max vertex for this simplex
		unsigned mid = (high + low)/2;
		if(bin.binom(mid, k) <= simplexHash) low = mid + 1;
		else high = mid;
	}
	return high - 1;
}

std::set<unsigned> simplexArrayList::getVertices(long long simplexHash, int dim, unsigned n){
	std::set<unsigned> v;
	for(unsigned k = dim+1; k>0; k--){ //Get all vertices by repeated binary search for max vertex
		n = maxVertex(simplexHash, n, k-1, k);
		v.insert(n);
		simplexHash -= bin.binom(n, k);
	}
	return v;
}

void simplexArrayList::prepareCofacets(int dim){
	indexConverter.clear();
	for(auto simplex : simplexList[dim+1]){
		indexConverter.insert(std::make_pair(simplex->hash, simplex));
	}
}

void simplexArrayList::prepareFacets(int dim){
	indexConverter.clear();
	for(auto simplex : simplexList[dim-1]){
		indexConverter.insert(std::make_pair(simplex->hash, simplex));
	}
}

std::vector<simplexNode_P> simplexArrayList::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<simplexNode_P, simplexNode_P>& pivotPairs, bool checkEmergent){
	std::vector<simplexNode_P> ret;
	int nPts = simplexList[0].size();
	unsigned k = simplex.size() + 1;
	std::set<unsigned>::reverse_iterator it = simplex.rbegin();
	long long index = simplexHash(simplex);

	//Try inserting other vertices into the simplex
	for(unsigned i=nPts; i-- != 0; ){
		if(it != simplex.rend() && i == *it - simplexOffset){ //Vertex i is already in the simplex
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

std::vector<simplexNode*> simplexArrayList::getAllCofacets(simplexNode_P simp, const std::unordered_map<long long, simplexNode_P>& pivotPairs, bool checkEmergent, bool recordVertices, unsigned dim){
	//Method builds out cofacets for incrementalPersistence

	std::vector<simplexNode*> ret;
	std::set<unsigned> vertices;

	if(recordVertices) vertices = simp->simplex;
	else vertices = getVertices(simp->hash, dim, simplexList[0].size());

	unsigned k = vertices.size() + 1;
	std::set<unsigned>::reverse_iterator it = vertices.rbegin();
	long long index = simp->hash;

	// Try inserting other vertices into the simplex
	for(unsigned i = simplexList[0].size(); i-- != 0; ){
		if(it != vertices.rend() && i == *it - simplexOffset){ //Vertex i is already in the simplex
			//Now adding vertices less than i -> i is now the kth largest vertex in the simplex instead of the (k-1)th
			index -= bin.binom(i, k-1);
			index += bin.binom(i, k); //Recompute the index accordingly
			--k;	//Now need to check for the (k-1)th vertex in the simplex
			++it; 	//Check for the previous vertex in the simplex (it is a reverse iterator)
		} else{
			double maxWeight = simp->weight;
			for(auto pt : vertices){
				maxWeight = std::max(maxWeight, (*distMatrix)[std::min(i, pt)][std::max(i, pt)]);
			}

			if(maxWeight <= maxEpsilon){ //Valid simplex
				simplexNode* x = new simplexNode();
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

std::vector<simplexNode*> simplexArrayList::getAllCofacets(simplexNode_P simp){
	return getAllCofacets(simp, std::unordered_map<long long, simplexNode_P>(), false, true, 0);
}

std::vector<simplexNode_P> simplexArrayList::getAllDelaunayCofacets(simplexNode_P simp){

	std::vector<simplexNode_P> ret;
	unsigned dimension  = simp->simplex.size();
        for(auto simplex : simplexList[dimension]){
                std::vector<unsigned> :: iterator it;
								std::vector<unsigned> v(simplex->simplex.size());
								it = std::set_intersection(simp->simplex.begin(),simp->simplex.end(),simplex->simplex.begin(),simplex->simplex.end(),v.begin());
								v.resize(it-v.begin());
								if(v.size() == simp->simplex.size())
									ret.push_back(simplex);
				}
	return ret;

}
std::vector<simplexNode*> simplexArrayList::getAllFacets(simplexNode* simp, bool recordVertices, unsigned dim){
	std::vector<simplexNode*> ret;
	std::set<unsigned> vertices;

	if(recordVertices) vertices = simp->simplex;
	else vertices = getVertices(simp->hash, dim + 1, simplexList[0].size());

	long long index = simp->hash;
	unsigned k = vertices.size();

	for(auto it = vertices.rbegin(); it != vertices.rend(); ++it){
		unsigned pt = *it;
		double maxWeight = 0;

		for(auto it = vertices.begin(); it != vertices.end(); ++it){
			if(*it == pt) continue;
			for(auto it2 = it; ++it2 != vertices.end(); ){
				if(*it2 == pt) continue;
				maxWeight = std::max(maxWeight, (*distMatrix)[*it][*it2]);
			}
		}

		simplexNode* x = new simplexNode();
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

std::vector<simplexNode*> simplexArrayList::getAllFacets(simplexNode_P simp, bool recordVertices, unsigned dimension){
	return getAllFacets(simp.get(), recordVertices, dimension);
}

std::vector<simplexNode_P> simplexArrayList::getAllFacets_P(simplexNode_P simp){
	std::vector<simplexNode_P> ret;

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
//
void simplexArrayList::insert(){
	//If this is the first point inserted...
	if(simplexList.size() == 0) simplexList.push_back({});

	unsigned i = simplexList[0].size();

	simplexNode_P insNode = std::make_shared<simplexNode>(simplexNode({i}, 0.0));
	insNode->hash = i;
	simplexList[0].insert(insNode);
}

// Search function to find a specific vector in the simplexArrayList
// weightedGraph[d][v][p] dimension d stores vectors v of point elements p of simplexes formed
bool simplexArrayList::find(std::set<unsigned> vector){
	if (simplexList.size() == 0) return false;

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

void simplexArrayList::initBinom(){
	bin = binomialTable(simplexList[0].size(), maxDimension+1);
}

// Expand the simplexArrayList to incorporate higher-level simplices
//	-> O(dnk) -> where n is the number of points, d is the dimension, and k is the number of d-1 simplices
//
//	Do this by comparing each simplex with points to insert
//
void simplexArrayList::expandDimensions(int dim){
	initBinom();

	//Iterate up to max dimension of simplex, starting at dim 2 (edges)
	for(unsigned d = 1; d <= dim; d++){

		//Check if we need to break from expanding dimensions (no more edges)
		if(simplexList.size() < d) break;
		if(simplexList.size() == d) simplexList.push_back({});

		//Iterate through each element in the current dimension's edges
		for(auto it = simplexList[d-1].begin(); it != simplexList[d-1].end(); it++){
			//Iterate over points to possibly add to the simplex
			//Use points larger than the maximal vertex in the simplex to prevent double counting
			unsigned minPt = *(*it)->simplex.rbegin() + 1;

			for(unsigned pt = minPt; pt < simplexList[0].size(); pt++){
				//Compute the weight using all edges
				double maxWeight = (*it)->weight;
				for(auto i : (*it)->simplex) maxWeight = std::max(maxWeight, (*distMatrix)[i][pt]);

				if(maxWeight <= maxEpsilon){ //Valid simplex
					simplexNode_P tot = std::make_shared<simplexNode>(simplexNode((*it)->simplex, maxWeight));
					tot->simplex.insert(pt);
					tot->hash = (*it)->hash + bin.binom(pt, tot->simplex.size());
					simplexList[d].insert(tot);
				}
			}
		}
	}
}

std::vector<simplexNode_P> simplexArrayList::expandDimension(std::vector<simplexNode_P> edges, bool recordVertices, unsigned dim){
	std::vector<simplexNode_P> nextEdges;

	//Iterate through each element in the current dimension's edges
	for(auto it = edges.begin(); it != edges.end(); it++){
		std::set<unsigned> vertices;

		if(recordVertices) vertices = (*it)->simplex;
		else vertices = getVertices((*it)->hash, dim - 1, simplexList[0].size());

		//Iterate over points to possibly add to the simplex
		//Use points larger than the maximal vertex in the simplex to prevent double counting
		unsigned minPt = *vertices.rbegin() + 1;

		for(unsigned pt = minPt; pt < simplexList[0].size(); pt++){
			//Compute the weight using all edges
			double maxWeight = (*it)->weight;
			for(auto i : vertices) maxWeight = std::max(maxWeight, (*distMatrix)[i][pt]);

			if(maxWeight <= maxEpsilon){ //Valid simplex
				simplexNode_P tot = std::make_shared<simplexNode>(simplexNode());
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

	if(recordVertices) std::sort(nextEdges.begin(), nextEdges.end(), cmpByWeight());
	return nextEdges;
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
//std::vector<double> :; simplexArrayList:: simplexHyperplane(std::set<unsigned> simplex,std::vector<std::vector<double>> inputData){

//}

void simplexArrayList:: buildAlphaComplex(std::vector<std::vector<int>> dsimplexmesh, int npts,std::vector<std::vector<double>> inputData){
unsigned maxDimension = dsimplexmesh[0].size()-1;
bin = binomialTable(npts,maxDimension+1);
for(int i=0;i<=maxDimension;i++)
	simplexList.push_back({});

for(auto simplex : dsimplexmesh){
	unsigned int pow_set_size = pow(2, simplex.size());
	for(int counter =1;counter<pow_set_size;counter++){
		double weight =0;
		std::set<unsigned> gensimp;
		for(int j=0;j<simplex.size();j++){
			if(counter & (1<<j)){
				unsigned indnew;
				indnew = *(std::next(simplex.begin(),j));
				for(auto x:gensimp){
					if(x<indnew){
						if(weight<(*distMatrix)[x][indnew])
							weight = (*distMatrix)[x][indnew];
					}
					else if(weight<(*distMatrix)[indnew][x])
						weight = (*distMatrix)[indnew][x];
				}
				gensimp.insert(indnew);
			}
		}
		//overwriting for Alpha
		double weight1;
		if(gensimp.size()>1)
	      	weight1 = utils::circumRadius(gensimp,distMatrix);
    else
		      weight1 = weight/2;

		simplexNode_P tot = std::make_shared<simplexNode>(simplexNode(gensimp,weight1));
		if(gensimp.size()>1)
				tot->circumRadius = utils::circumRadius(gensimp,distMatrix);
		else{
		    tot->circumRadius = weight/2;
			}
		if(gensimp.size()>2)
				tot->circumCenter = utils::circumCenter(gensimp,inputData);
	//	tot->simplexhyperplane = simplexHyperplane(gensimp,inputData);

		else if(gensimp.size()==2){
 			auto first = gensimp.begin();
			std::vector<double> R;
			std::vector<double> A = inputData[*first];
      std::advance(first, 1);
			std::vector<double> B = inputData[*first];
   		std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
		  tot->circumCenter = R;
    }
		else
   		tot->circumCenter = inputData[*(gensimp.begin())];


		if(gensimp.size()==1)
			tot->hash = *(gensimp.begin());
		else
			tot->hash = simplexHash(gensimp);
		simplexList[gensimp.size()-1].insert(tot);
		gensimp.clear();
	}
}

//ALPHA COMPLEX FILTERTION BASED on FOLLOWING algorithm
/*Filtration value computation algorithm


for i : dimension →0 do
   for all σ of dimension i
        if filtration(σ) is NaN then
            filtration(σ)=α2(σ)
        end if
        for all τ face of σ do            // propagate alpha filtration value
          if  filtration(τ) is not NaN then
               filtration(τ) = min( filtration(τ), filtration(σ) )
          else
             if τ is not Gabriel for σ then
               filtration(τ) = filtration(σ)
          end if
       end if
     end for
  end for
end for

make_filtration_non_decreasing()
prune_above_filtration()

  */

for(int dim = simplexList.size()-1;dim>=0;dim--){
	   for(auto simplexiter = simplexList[dim].rbegin();simplexiter != simplexList[dim].rend();simplexiter++){
			    simplexNode_P simplex = (*simplexiter);
	        if(simplex->filterationvalue ==-1)
	            simplex->filterationvalue = simplex->circumRadius;
				//	for(int facedim = dim;facedim>=0;facedim--)
				      if(dim>0){
							for(auto face :simplexList[dim-1]){
								bool gabriel = true;
								std::vector<unsigned> points_check(simplex->simplex.size());
								std::vector<unsigned> guilty_points_check;
								std::vector<unsigned> :: iterator it;
								std::vector<unsigned> v(face->simplex.size());
								it = std::set_intersection(simplex->simplex.begin(),simplex->simplex.end(),face->simplex.begin(),face->simplex.end(),v.begin());
								v.resize(it-v.begin());
								if(v.size() == face->simplex.size()){
                     if(face->filterationvalue !=-1){
										 			face->filterationvalue = std::min(face->filterationvalue ,simplex->filterationvalue);
												}
										 else {
													std::vector<unsigned>::iterator it;
													it=std::set_difference (simplex->simplex.begin(), simplex->simplex.end(), face->simplex.begin(), face->simplex.end(), points_check.begin());
													points_check.resize(it-points_check.begin());
    											for(it = points_check.begin(); it !=points_check.end();++it){
														std::vector<double> coordinates;
                            for(int i =0;i<face->circumCenter.size();i++)
																	coordinates.push_back(inputData[*it][i]);
														double distance = utils::vectors_distance(coordinates,face->circumCenter);
														if(pow(distance,2)<face->circumRadius){
															 gabriel = false;
															 guilty_points_check.push_back((*it));
													}
		   							}
						  		}

							}
							if(!gabriel){
										std::vector<unsigned> v(simplex->simplex.size());
										std::vector<unsigned>::iterator it;
										it=std::set_union (face->simplex.begin(), face->simplex.end(), guilty_points_check.begin(), guilty_points_check.end(), v.begin());
										v.resize(it-v.begin());
										for(auto face1 :simplexList[v.size()-1]){
											std::vector<unsigned> v1(simplex->simplex.size());
											it = std::set_intersection(v.begin(),v.end(),face1->simplex.begin(),face1->simplex.end(),v1.begin());
											v1.resize(it-v1.begin());
											if(v1.size() == face1->simplex.size())
														face->filterationvalue = face1->filterationvalue;

										   }
					           }
									}
					}
   }
}
//Reinserting to sort by filterationvalue and remove simplexes with weight greater than alphafilteration value (maxEpsilon)
std::vector<std::set<simplexNode_P, cmpByWeight>> simplexList1;		//Holds ordered list of simplices in each dimension
for(int dim=0;dim < simplexList.size();dim++){
	   simplexList1.push_back({});
	   for(auto simplex :simplexList[dim]){
			 if(simplex->filterationvalue <= maxEpsilon){ //Valid Simplex after filteration
			 			 simplex->weight = simplex->filterationvalue;
			 		   simplexList1[dim].insert(simplex);
					 }
		 }
	 }
simplexList = simplexList1;
return;
}



std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> simplexArrayList::recurseReduce(simplexNode_P simplex, std::vector<std::set<unsigned>> removals, std::vector<std::set<unsigned>> checked){
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

simplexArrayList::~simplexArrayList(){
	simplexList.clear();
}
