#pragma once

// Header file for bettiPipe class - see bettiPipe.cpp for descriptions
#include <map>
#include <vector>
#include "basePipe.hpp"
#include "simplexBase.hpp"
#include "simplexArrayList.hpp"
#include "utils.hpp"

class involutedPersistence : public basePipe {
	private:
		std::vector<simplexNode*> edges;
		unsigned dim;

		struct sortLexicographic{ //Sort nodes by weight, then by lexicographic order
			bool operator()(simplexNode* a, simplexNode* b) const{
				if(a->weight == b->weight){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
					auto itA = a->simplex.begin(), itB = b->simplex.begin();
					while(itA != a->simplex.end()){
						if(*itA != *itB) return *itA < *itB;
						++itA; ++itB;
					}
					return false;
				} else{
					return a->weight < b->weight;
				}
			}
		};

	public:
	    involutedPersistence();
	    void setupSimplices(std::vector<simplexNode*>, unsigned);
	    void runPipe(pipePacket &inData);
	    bool configPipe(std::map<std::string, std::string> &configMap);
		void outputData(pipePacket&);

		template <class simplexNodePointer, class comp>
		std::vector<simplexNodePointer> test(pipePacket& inData, std::vector<simplexNodePointer> edges, std::vector<simplexNodePointer> pivots, unsigned dimension, comp compStruct, std::string mode){
			std::sort(edges.begin(), edges.end(), compStruct());

			std::sort(pivots.begin(), pivots.end(), compStruct);
			typename std::vector<simplexNodePointer>::iterator it = pivots.begin();


			std::vector<simplexNode_P> nextPivots;	 					//Pivots for the next dimension
			std::unordered_map<simplexNodePointer, std::vector<simplexNodePointer>> v;				//Store only the reduction matrix V and compute R implicity
			std::unordered_map<long long, simplexNodePointer*> pivotPairs;						//For each pivot, which column has that pivot

			//Iterate over columns to reduce in reverse order
			for(auto columnIndexIter = edges.begin(); columnIndexIter != edges.end(); columnIndexIter++){
				simplexNodePointer simplex = (*columnIndexIter);		//The current simplex

				//Only need to test hash for equality
				//Not a pivot -> need to reduce
				if((*it)->hash != simplex->hash){

					//Get all cofacets using emergent pair optimization
					std::vector<simplexNode*> faceList = (mode == "homology" ? inData.complex->getAllFacets(simplex) : inData.complex->getAllCofacets(simplex));

					std::vector<simplexNodePointer> columnV;	//Reduction column of matrix V
					columnV.push_back(simplex); //Initially V=I -> 1's along diagonal

					//Build a heap using the coface list to reduce and store in V
					std::make_heap(faceList.begin(), faceList.end(), compStruct());

					while(true){
						simplexNode* pivot;

						while(!faceList.empty()){
							pivot = faceList.front();

							//Rotate the heap
							std::pop_heap(faceList.begin(), faceList.end(), compStruct());
							faceList.pop_back();

							if(!faceList.empty() && pivot->hash == faceList.front()->hash){ //Coface is in twice -> evaluates to 0 mod 2
								if(inData.complex->simplexType == "simplexArrayList"){
									delete pivot;
									delete faceList.front();
								}

								//Rotate the heap
								std::pop_heap(faceList.begin(), faceList.end(), compStruct());
								faceList.pop_back();
							} else{

								faceList.push_back(pivot);
								std::push_heap(faceList.begin(), faceList.end(), compStruct());
								break;
							}
						}

						if(faceList.empty()){ //Column completely reduced
							break;
						} else if(pivotPairs.find(pivot->hash) == pivotPairs.end()){ //Column cannot be reduced
							pivotPairs.insert({pivot->hash, simplex});
							nextPivots.push_back(std::shared_ptr<simplexNode>(pivot));

							std::sort(columnV.begin(), columnV.end());
							auto it = columnV.begin();
							while(it != columnV.end()){
								if((it+1) != columnV.end() && *it==*(it+1)) ++it;
								else v[simplex].push_back(*it);
								++it;
							}

							if(simplex->weight != pivot->weight){
								bettiBoundaryTableEntry des = { dim, pivot->weight, simplex->weight, ut.extractBoundaryPoints(v[simplex]) };
								inData.bettiTable.push_back(des);
							}

							//Don't delete the first entry because that is converted to a smart pointer and stored as a pivot
							if(inData.complex->simplexType == "simplexArrayList"){
								for(int i=1; i<faceList.size(); i++) delete faceList[i];
							}

							break;
						} else{ //Reduce the column of R by computing the appropriate columns of D by enumerating cofacets
							for(simplexNodePointer simp : v[pivotPairs[pivot->hash]]){
								columnV.push_back(simp);
								std::vector<simplexNode*> faces = (mode == "homology" ? inData.complex->getAllFacets(simp) : ((simplexArrayList*) inData.complex)->getAllCofacets(simp));
								faceList.insert(faceList.end(), faces.begin(), faces.end());
							}
							std::make_heap(faceList.begin(), faceList.end(), sortLexicographic());
						}
					}
				} else ++it;
			}

			return nextPivots;
		}
};

