#include <string>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <iostream>
#include <map>
#include "indSimplexTree.hpp"


indSimplexTree::indSimplexTree(double _maxEpsilon, std::vector<std::vector<double>> _distMatrix, int _maxDim){
	indexCounter = 0;
	distMatrix = _distMatrix;
	maxDim = _maxDim;
	maxEpsilon = _maxEpsilon;
	simplexType = "indSimplexTree";
	return;
}

std::pair<std::vector<std::set<unsigned>>,std::vector<std::set<unsigned>>> indSimplexTree::recurseReduce(std::set<unsigned> sourceSet, int fc, std::set<unsigned> curNode,int d, std::vector<std::set<unsigned>> removalSimplices, std::vector<std::set<unsigned>> processedSimplices){
	processedSimplices.push_back(curNode);
	
	std::vector<std::set<unsigned>> localBuf;
	std::vector<std::set<unsigned>> localFaces;
	
	//Get simplex subsets to track common faces - these are candidate subsets before evaluation
	auto subsets = ut.getSubsets(curNode, d);
	
	//For each candidate subset prior to evaluation
	for(auto subset : subsets){
		for(int z = 0; z < dimensions[d].size(); z++){
			
			//Check if the simplex has already been processed
			if(std::find(removalSimplices.begin(), removalSimplices.end(), dimensions[d][z]->simplexSet) == removalSimplices.end() && \
				std::find(processedSimplices.begin(), processedSimplices.end(), dimensions[d][z]->simplexSet) == processedSimplices.end()){
			//if(dimensions[d][z]->simplexSet != curNode){
				
				
				//If there is an intersection with the current simplex
				if(ut.setIntersect(dimensions[d][z]->simplexSet,subset,true).size() == dimensions[d][z]->simplexSet.size() - 1	){
					localBuf.push_back(subset);
					localFaces.push_back(dimensions[d][z]->simplexSet);
					//Add the current simplex to graph removal
					//removalSimplices.push_back(dimensions[d][z]->simplexSet);
					//std::cout << "Removal: ";
					//ut.print1DVector(dimensions[d][z]->simplexSet);
					//deletion(dimensions[d][z]->simplexSet);
									
					
					//Find the additional simplices to remove from the graph					
					//if(std::find(removalSimplices.begin(), removalSimplices.end(), curNode) == removalSimplices.end()){
					//	removalSimplices.push_back(curNode);
					//}
					break;
				}
			}
		}
	}
	//std::cout << "Local Buf: " << localBuf.size() << "\tSubsets: " << subsets.size() << std::endl;
	
	//Recurses to determine if the neighboring faces all have connected faces
	if(localBuf.size() == subsets.size()){
		bool recurse = false;
		for(auto face : localFaces){
			
			auto a = recurseReduce(curNode, 0, face, d, removalSimplices, processedSimplices);
			processedSimplices = a.first;
			if(a.second.size() > removalSimplices.size())
				recurse = true;
			removalSimplices = a.second;
			
		}
		
		if(recurse){
			//Refresh the local buffer for faces that match
			auto a = recurseReduce(sourceSet, 0, curNode, d, removalSimplices, processedSimplices);
		}
		
		
	} else if(dimensions[d].size() == 1 || localBuf.size() != subsets.size()){		
		
		std::set<unsigned> non_min;
		double minVal = -1;
		
		//For each subset
		for(auto z : subsets){
			//if the subset is not in the local buffer (i.e. did not get filtered by faces)
			if(std::find(localBuf.begin(), localBuf.end(), z) == localBuf.end()){
				//std::cout << "Candidate: ";
				//ut.print1DVector(z);
				double d = getWeight(z);
				if(minVal < 0){
					minVal = d;
					non_min = z;
				} else if(d > minVal){
					non_min = z;
				} else {
					minVal = d;
				}
				//std::cout << "\tWT: " << d << std::endl;
			}
		}		
		if(non_min.size() > 0){
			removalSimplices.push_back(non_min);
			removalSimplices.push_back(curNode);
			std::cout << "Removal: ";
			ut.print1DVector(non_min);
			//std::cout << "  ,  ";
			ut.print1DVector(curNode);
			//auto a = recurseReduce(curNode, 0, non_min, d-1, removalSimplices, processedSimplices);
			//processedSimplices = a.first;
			//removalSimplices = a.second;
		}
		fc = 0;
	}
	
	return std::make_pair(processedSimplices, removalSimplices);
}

std::vector<std::vector<simplexBase::graphEntry>> indSimplexTree::coreduction(graphEntry cur){
	std::vector<std::vector<graphEntry>> ret = indexedGraph;
	//Queue q ; q.insert (cur)
	std::vector<graphEntry> simplexQueue;
	
	simplexQueue.push_back(cur);
		
	//While (queue.next() != end)
	while(simplexQueue.size() > 0){
		
			//std::cout << "cored" << std::endl;
		// s = q.pop();
		graphEntry currentEntry = simplexQueue[0];
		
		//Compute cbds, bds
		//	cbds(a) = {s <= S | k(s,t) != 0 for some t <= a}
		std::vector<graphEntry> cbds;
		
		//	bds(a) = {t <= S | k(s, t) != 0 for some s <= a}
		std::vector<graphEntry> bds;
		
		
			//std::cout << "cored " << simplexQueue.size() << std::endl;
		//if bds.size = 1
		if(bds.size() == 1){
		
			//remove s from indexedGraph
			
			
			//foreach coboundary of t
			for(auto ge : cbds){
				//queue.insert(u)
				simplexQueue.push_back(ge);
			}
			
			//remove t from s
			
		}
		
		
		//if bds s = null 
		if(bds.size() == 0){
			
			//for each coboundary of s
			for(auto ge : cbds){
				//queue.insert(u)
				simplexQueue.push_back(ge);
			}
		}
		
		simplexQueue.erase(simplexQueue.begin());
		//std::cout << "cored " << simplexQueue.size() << std::endl;
	}
	//std::cout << "RET" << std::endl;
	return ret;
}

void indSimplexTree::recurseInsert(indTreeNode* node, unsigned curIndex, int depth, double maxE, std::set<unsigned> simp){
	//Incremental insertion
	//Recurse to each child (which we'll use the parent pointer for...)
	indTreeNode* temp;
	
	double curE = distMatrix[node->index][indexCounter];
	curE = curE > maxE ? curE : maxE;
	
	depth = simp.size() - 1;
	
	
	//Check if the node needs inserted at this level
	if(curE < maxEpsilon){
		indTreeNode* insNode = new indTreeNode();
		insNode->index = curIndex;
		nodeCount++;
		std::set<unsigned> newSimp = simp;
		newSimp.insert(node->index);
		insNode->simplexSet = newSimp;
		dimCounts[depth+1]++;
		
		//Get the largest weight of this simplex
		maxE = curE > node->weight ? curE : node->weight;
		insNode->weight = maxE;
		
		//Check the dimensional list
		//	if no nodes exist at the dimension, this is the first
		if(depth == dimensions.size() - 1){
			std::vector<indTreeNode*> tn;
			tn.push_back(insNode);
			dimensions.push_back(tn);
		} else {
			dimensions[depth+1].push_back(insNode);
		}
		
		
		//Check if the node has children already... (remember parent == child for incremental)
		if(node->child == nullptr){
			insNode->parent = node;
			node->child = insNode;
		} else {
			insNode->parent = node;
			insNode->sibling = node->child;
			node->child = insNode;
			
			temp = insNode->sibling;
			//Have to check the children now...
			if(newSimp.size() <= maxDim){
				do {
					recurseInsert(temp, curIndex, depth + 1, maxE, newSimp);
				} while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);
			}
		}
	}
		
	return;
}

void indSimplexTree::printTree(indTreeNode* head){
	
	std::cout << "_____________________________________" << std::endl;
	
	if(head == nullptr){
		std::cout << "Empty tree... " << std::endl;
		return;
	}
	
	indTreeNode* current;
	
	for(int i = 0; i < dimensions.size() ; i++){
		std::cout << std::endl << "Dim: " << i << std::endl << std::endl;
		
		current = dimensions[i][0];
		
		do{
			
			std::cout << current->index << "," << current << "," << current->parent << "\t";
			
		} while(current->sibling != nullptr && (current = current->sibling) != nullptr);
		
		std::cout << std::endl;
	}
	
	std::cout << "_____________________________________" << std::endl;
	return;
}
	
	
void indSimplexTree::insertInductive(){
	//Create our new node to insert
	indTreeNode* curNode = new indTreeNode;
	curNode->index = indexCounter;
	
	//Loop each dimensional list; start at dn, then dn-1, ..., d0 (OLD WAY, INDUCTIVE VR)
	for(int i = (dimensions.size()-1 < maxDim ? dimensions.size() - 1 : maxDim); i >= 0; i--){
		//Loop each sibling in the current dimensional list (as a do-while, break when no more siblings)
		curNode = dimensions[i][0];
		do{
			//Determine if the index needs to inserted below the current node
			//	First, check the distance matrix for curNodeIndex v. curIndex
			//	Second (if first is true), check all parents dist matrix of parentNodeIndex v. curIndex
			if(distMatrix[curNode->index][indexCounter] < maxEpsilon){
				//	This node is a candidate, now we need to check each parent of the current branch
				//		to ensure the distance matrix entry for parentNodeindex v. curIndex is < maxEpsilon

				bool ins = true;
				indTreeNode* parentNode = curNode; //Temporarily store the parent so we don't lose track of our current
			
				// This loop will check each parent and set whether to insert a node
				do{
					if(distMatrix[parentNode->index][indexCounter] > maxEpsilon){
						ins = false;
					}	
				}while(ins && parentNode->parent != nullptr && (parentNode = parentNode->parent) != nullptr);
				
				//Insert the node; the distance to each parent is less than epsilon
				if(ins){
					//Allocate a new node to be inserted into the tree
					indTreeNode* insNode = new indTreeNode;
					
					//This new node's parent is the current node we're indexing on
					insNode->parent = curNode;
					insNode->index = indexCounter;
					
					insNode->sibling = nullptr;
				
					//Check the dimensional list
					//	if no nodes exist at the dimension, this is the first
					if(i == dimensions.size() - 1){
						
						std::vector<indTreeNode*> tempInd;
						tempInd.push_back(insNode);
						dimensions.push_back(tempInd);				
					
					//	Nodes currently exist at the dimension, insert as a sibling
					//		This should be inserted with nodes of the same parent
					} else {
						indTreeNode* iterateNode = dimensions[i+1][0];
						ins = false; //Store if we've inserted this node while looping through siblings
						
						//Loop through dimensional node siblings
						do{
								//Check if we've found parent nodes...
								if(iterateNode->parent == insNode->parent){
									ins = true;
									
									//perform an insertion of the node...
									insNode->sibling = iterateNode->sibling;
									iterateNode->sibling = insNode;
									
									nodeCount++;
								}
						} while(!ins && iterateNode->sibling != nullptr && (iterateNode = iterateNode->sibling) != nullptr);
							
						// Check if we successfully inserted a node; if not, insert to end
						if(!ins){
							iterateNode->sibling = insNode;
						}
					}
				}
			}			
			
		}while(curNode->sibling != nullptr && (curNode = curNode->sibling) != nullptr);
	}
	
	return;
}




// Insert a node into the tree
//		
void indSimplexTree::insert(std::vector<double>&) {
	
	if(distMatrix.size() == 0){
		ut.writeDebug("simplexTree","Distance matrix is empty, skipping insertion");
		return;
	}
	
	//Create our new node to insert
	indTreeNode* curNode = new indTreeNode;
	curNode->index = indexCounter;
	
	std::set<unsigned> tempSet = {curNode->index};
	
	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(head == nullptr){	
		
		curNode->simplexSet = { indexCounter };	
		curNode->sortedIndex = indexCounter;
		head = curNode;
		indexCounter++;
		nodeCount++;
		std::vector<indTreeNode*> tempInd;
		tempInd.push_back(head);
		dimensions.push_back(tempInd);
		
		dimCounts[0]++;
		
		//std::vector<indGraphEntry> tempWEG;
		//tempWEG.push_back(new graphEntry(tempSet, 0, curNode));
		//indexedGraph.push_back(tempWEG);
		
		return;
	}
		
		
	// This needs to be a recursive span -> (or not!)
	//		if the node has a child, recurse
	//		if the node has a sibling, recurse
	//
	//d0 -->          ({!})
	//			       /
	//			      /
	//d1 -->	   | 0 | 1 | ... |
	//		        /
	//             /
	//d2 -->    | 1 | 2 | 3 | 
	//           /         \
	//	        /           \
	//d3 --> | 2 | 3 |     | 4 | 5 |
	//

	indTreeNode* temp = dimensions[0][0];
	
	
	
	//Now let's do it the new way - incremental VR;
	//	Start at the first dimension (d0);
	//		if e (d0, dn) < eps
	//			recurse down tree
	//				recurse	
	//			insert to current
	//	iterate to d0->sibling
	
	
	
	do{
		recurseInsert(temp, indexCounter, 0, 0, tempSet);
	}while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);
	
	//Insert into the right of the tree
	temp = dimensions[0][0];
	indTreeNode* ins = new indTreeNode();
	while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);
	ins->index = indexCounter;
	ins->sibling = nullptr;
	ins->parent = nullptr;
	ins->sortedIndex = indexCounter;
	ins->simplexSet = { indexCounter };
	temp->sibling = ins;
	//indexedGraph[0].push_back(std::make_pair(ins, std::make_pair(tempSet, 0)));
	
	
	dimensions[0].push_back(ins);
	dimCounts[0]++;
	nodeCount++;
	indexCounter++;
	return;

}


// Iterative function to search for a key in the tree. The function returns true
// if the key is found, else it returns false.
bool indSimplexTree::search(std::set<unsigned> simplex){
	indTreeNode* curNode = dimensions[0][0];
	
	for(auto i :simplex){
		while(curNode != nullptr && curNode->index != i);
		if(curNode != nullptr){ curNode = curNode->child; }
		else { return false; }
	}
	return true;
}


double indSimplexTree::getWeight(std::set<unsigned> search){
	indTreeNode* curNode = dimensions[0][0];
	
	//std::cout << "Get Weight - Size: " << search.size() << "\tSimplex Count: " << dimensions[search.size()-1].size() << "\t";
	//ut.print1DVector(search);
	
	for(auto a : dimensions[search.size()-1]){
		if(search == a->simplexSet){
			return a->weight;
		}
	}
	return 10000;
	
	/*	
	for(auto i : dimensions[search.size()]){
		if
		while(curNode != nullptr && curNode->index != i){curNode = curNode->sibling;};
		if(curNode != nullptr){curNode = curNode->child;}
		else return 10000;
	}
	if(curNode != nullptr){
		ut.print1DVector(curNode->simplexSet);
		std::cout << "node " << curNode->weight << std::endl;
		return curNode->weight;
	} else {
		return 10000;
	}*/
}


// The following function returns true if a given node has a child.
bool indSimplexTree::haveChild(indSimplexTree const* present) {

}



// A recursive function to delete a simplex (and sub-branches) from the tree.
bool indSimplexTree::deletion(std::set<unsigned> removalEntry) {
	//Remove the entry in the simplex tree
	bool found = true;
	indTreeNode* curNode = dimensions[0][0];
	int d = removalEntry.size() - 1;
	
	
	for(int i = 0; i < dimensions[d].size(); i++){
		if(removalEntry == dimensions[d][i]->simplexSet){
			curNode = dimensions[d][i];
			dimensions[d].erase(dimensions[d].begin() + i);
			return deletion(curNode);
		}
	}
	return false;
			
			
			/*
	auto subsets = ut.getSubsets(curNode, d);
	for(int i = 0; i < subsets.size(); i++){
		if(subsets[i] == sourceSet)
			subsets.erase(subsets.begin() + i);
		
		//Find the current entry by index
		while(curNode!= nullptr && curNode->index != i){curNode = curNode->sibling;};
		
		
		if(curNode!=nullptr){curNode = curNode->child;}
		
		else{
			found = false;
			break;
		}
	}
	
	std::cout << "DELETION!" << std::endl;
	return ;
	
	*/
	
	
}


// A recursive function to delete a simplex (and sub-branches) from the tree.
bool indSimplexTree::deletion(indTreeNode* removalEntry) {
	indTreeNode* curNode = removalEntry;
	
	//Iterate to the bottom of branch
	//while(curNode != nullptr){ 
	//	curNode->child->parent = curNode;
	//	curNode=curNode->child;
	//}
	
	
	//If we did go down, remove on the way back up
	while(curNode != removalEntry){ 
		curNode = curNode->parent;
		delete curNode->child;
		curNode->child = nullptr;
	}
	
	//curNode = curNode->parent;
	delete curNode;
	nodeCount--;
	//curNode->child = nullptr;
	
	//NOTE: Also need to remove from the indexed graph (if done after shuffle)
	
	
	return false;
}

int indSimplexTree::vertexCount(){
	//Return the number of vertices in the tree
	return nodeCount;
}

int indSimplexTree::simplexCount(){
	//Return the number of simplices in the tree
	return nodeCount;
} 

bool indSimplexTree::find(std::set<unsigned> simplex){
	if(dimensions.size() == 0){
		ut.writeDebug("indSimplexTree","Complex is empty, skipping find");
		return false;
	}
		
	indTreeNode* curNode = dimensions[0][*simplex.begin()];
	
	for(auto i = simplex.begin() ; i != simplex.end(); i++){
		if(i != simplex.begin() && curNode->child != nullptr){
			curNode = curNode->child;
			
			while(curNode != nullptr){
				if(curNode->index != *i){
					curNode = curNode->sibling; 
				} else if(curNode->index == *i){
					break;
				} else { 
					return -1; 
				}
			}
		}
	}
	return curNode->sortedIndex;	
	
}

double indSimplexTree::getSize(){
	//Size of node: [int + byte (*) + byte (*)] = 18 Bytes
	return nodeCount * (sizeof(indTreeNode) + sizeof(graphEntry));
}

bool indSimplexTree::compareByWeight(const graphEntry &a, const graphEntry &b){
	return a.weight < b.weight;	
}

void indSimplexTree::sortAndBuildGraph(){
	
	std::cout << "========================    SORTING    ===============================" << std::endl;
	indTreeNode* curNode = dimensions[0][0];
	std::vector<graphEntry> curEntry;
	std::vector<std::vector<graphEntry>> ret;
	
	//Handle points first, default ordering
	int i = 0;
	for(auto cur : dimensions[0]){
		graphEntry ge(cur->simplexSet, cur->weight, cur);
		curEntry.push_back(ge);
		//ut.print1DVector(cur->simplexSet);
		i++;
	}	
	ret.push_back(curEntry);
	//std::cout << "Sort and Build " << dimensions.size() << std::endl;
	
	//When inserting the next few dimensions, need to first update 
	//	the indTreeNode to point to the correct graphEntry now (since we've resorted)
	//
	//	This just involves iterating the sorted vector and following the indTreeNode
	//		pointer to update the curNode; we'll use this to speed up checkFace functions

	int iter = 0;
	for(int i = 1; i < dimensions.size(); i++){
		
		std::vector<graphEntry> innerEntry;//(dimCounts[i], *(new graphEntry()));
		
		for(auto cur : dimensions[i]){
			indSimplexTree::graphEntry ge(cur->simplexSet, cur->weight, cur);			
			innerEntry.push_back(ge);
		}		
		
		
		//SORT graph vector HERE, reassociate indTreeNode pointer
		std::sort(innerEntry.begin(), innerEntry.end(), indSimplexTree::compareByWeight);
		//std::cout << "Sorted!" << std::endl;
		
		int it = 0;
		for(auto a : innerEntry){
			a.entry->sortedIndex = it;
			it++;
			//std::cout << a.weight << "\t";
			//ut.print1DVector(a.simplexSet);
		}
		
		//std::cout << "Pushing..." << std::endl;
		
		ret.push_back(innerEntry);
	}
	
	indexedGraph = ret;
	
	
}


std::vector<std::vector<std::pair<std::set<unsigned>,double>>> indSimplexTree::getAllEdges(double epsilon){
	
	
	
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> ret;
	
	if(dimensions.size() == 0){
		ut.writeDebug("indSimplexTree","Complex is empty, no edges to return");
		return ret;
	}
	
	if(!isSorted){
		sortAndBuildGraph();
		isSorted = true;
	}
	
	for(std::vector<graphEntry> a : indexedGraph){
		std::vector<std::pair<std::set<unsigned>,double>> tempEdges;
		for(auto b : a){
			tempEdges.push_back(std::make_pair(b.simplexSet, b.weight));
		}
	
		ret.push_back(tempEdges);
	}
	
	return ret;
	
}



std::vector<std::vector<indSimplexTree::graphEntry>> indSimplexTree::getIndexEdges(double epsilon){
	std::vector<std::vector<indSimplexTree::graphEntry>> ret;
	
	if(dimensions.size() == 0){
			ut.writeDebug("indSimplexTree","Complex is empty, no edges to return");
			return ret;
	}

	if(!isSorted){
		sortAndBuildGraph();
		isSorted = true;
	}
	return indexedGraph;
}

void indSimplexTree::reduceComplex(){
	std::vector<graphEntry> curEntry;
	graphEntry ge;
	
	std::cout << "Reducing simplex tree elementary collapse" << std::endl;
	std::cout << "\tOriginal simplex count: " << nodeCount << std::endl;
	
	if(maxDim < 2){
		std::cout << "Dim < 2, exiting reduction" << std::endl;
		return;
	}
	
	
	//Start at the highest dimension and check for cofaces of simplices
	//		if there isn't a shared face with another simplex of the
	//		same dimension then remove the face
	//
	//
	//	There are 4 cases that can happen as we iterate a simplex s_d:
	//
	//		1. All faces are shared, s_d >= max_d
	//		2. All faces are shared, s_d < max_d
	//		3. Some faces are shared
	//		4. No faces shared
	//
	//
	//	To do this, we recurse through shared faces and examine neighbor simplices;
	//
	//		If a neighbor simplex has all shared faces (and recurses/faces the same), no removal simplices
	//			will be added. This will be case 1 / 2
	//
	//		If a neighbor simplex does not have shared faces, remove the largest edge and the simplex 
	//			itself. Continue iterating other faces of original simplex. This will be case 3.
	//
	//		After all faces are finished, if processed simplices.size() = 1, case 4.
	//
	
	//Iterate each dimension of simplices
	for(int d = dimensions.size()-1; d > 1 ; d--){
		std::vector<graphEntry> curGraph;
		std::map<std::set<unsigned>,int> counts;
		std::vector<std::set<unsigned>> removalSimplices;
		std::vector<std::set<unsigned>> processedSimplices;
		
		
		//Iterate each node in the simplex list, dimension d
		for(auto curNode : dimensions[d]){	
			
			//Check if the current node has been processed already
			if(std::find(processedSimplices.begin(), processedSimplices.end(), curNode->simplexSet) == processedSimplices.end()){
				
				//If the current node hasn't been processed, recursively iterate through shared faces and remove simplices
				//		This function handles case 3 and 4.
				auto a =  recurseReduce({}, 0, curNode->simplexSet, d, removalSimplices, processedSimplices);
				
				//Check for cases 1 / 2
				if(removalSimplices.size() == a.second.size()){
					
					// 1. All faces are shared, s_d >= max_d
					if(curNode->simplexSet.size() >= maxDim){
						// Do nothing, need to process through t-array if a feature forms 					

					// 2. All faces are shared, s_d < max_d
					} else {
						//Use the minimum path algorithm to remove edges from processed simplices
						//	Need to remove all non-minimal faces of each processed simplex, along
						//		with the simplex itself
					}
				}
				
				processedSimplices = a.first;
				removalSimplices = a.second;
			}
					
		}
		
		//std::cout << "_________Removals_________" << std::endl;
		//for(auto p: removalSimplices){
		//	std::cout << "\t";
		//	ut.print1DVector(p);
		//}
		//std::cout << std::endl;
		
		//Remove the removals...
		for(auto r : removalSimplices){
			for(int i = 0; i < dimensions[d].size(); i++){
				if(r == dimensions[d][i]->simplexSet){
					dimensions[d].erase(dimensions[d].begin() + i);
					nodeCount--;
				}
			}
			deletion(r);
		}
			
		
	}
	
	
	std::cout << "\tElementary collapse count: " << nodeCount << std::endl;
	std::cout << "Reducing simplex tree coreduction" << std::endl;
	
	
	
	for(auto a : indexedGraph){
		for(auto b: a){
			
			//	Elementary coreduction pair if bds(b) = {a}
			//
			//	bds(a) = {t <= S | k(s, t) != 0 for some s <= a}
			//		where dim(s) = dim(t) + 1
			//
			
			coreduction(b);
			//std::cout << "Reduced!" << std::endl;
			
		}
	}
	
	std::cout << "\tCoreduction count: " << nodeCount << std::endl;
		
	return;
}


std::vector<std::vector<unsigned>> indSimplexTree::getDimEdges(int,double){
	std::vector<std::vector<unsigned>> ret;
	return ret;
	
}
