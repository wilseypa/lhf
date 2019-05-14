#include <string>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <iostream>
#include "indSimplexTree.hpp"
#include "utils.hpp"


bool sortSecond(const std::pair<std::set<unsigned>, double> &a, const std::pair<std::set<unsigned>, double> &b){
	return (a.second < b.second);
}

indSimplexTree::indSimplexTree(double _maxEpsilon, std::vector<std::vector<double>> _distMatrix, int _maxDim){
	indexCounter = 0;
	distMatrix = _distMatrix;
	maxDim = _maxDim;
	maxEpsilon = _maxEpsilon;
	simplexType = "indSimplexTree";
	return;
}

void indSimplexTree::recurseInsert(indTreeNode* node, unsigned curIndex, int depth, double maxE, std::set<unsigned> simp){
	utils ut;
	//Incremental insertion
	//Recurse to each child (which we'll use the parent pointer for...)
	indTreeNode* temp;
	
	double curE = distMatrix[node->index][indexCounter];
	curE = curE > maxE ? curE : maxE;
	
	
	//Check if the node needs inserted at this level
	if(curE < maxEpsilon){
		indTreeNode* insNode = new indTreeNode();
		insNode->index = curIndex;
		nodeCount++;
		simp.insert(node->index);
		
		//Get the largest weight of this simplex
		maxE = curE > node->weight ? curE : node->weight;
		insNode->weight = maxE;
		
		//if depth (i.e. 1 for first iteration) is LT weightGraphSize (starts at 1)
		//std::cout << "WEG: " << weightEdgeGraph.size() << "\t" << simp.size() << "\t" << maxDim << std::endl;
		
		//graphEntry ge;
		
		//if(indexedGraph.size() < simp.size()){
		//	std::vector<graphEntry> tempWEG;
		//	ge.weight = maxE;
		//	ge.simplexSet = simp;
		//	tempWEG.push_back(ge);
		//	indexedGraph.push_back(tempWEG);
		//} else {
		//	indexedGraph[simp.size() - 1].push_back(std::make_pair(insNode, std::make_pair(simp, maxE)));
		//}
		
		//Check if the node has children already... (remember parent == child for incremental)
		if(node->parent == nullptr){
			node->parent = insNode;
			
		} else {
			insNode->sibling = node->parent;
			node->parent = insNode;
			
			temp = insNode->sibling;
			//Have to check the children now...
			if(simp.size() <= maxDim){
				do {
					recurseInsert(temp, curIndex, depth + 1, maxE, simp);
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
		
		current = dimensions[i];
		
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
		curNode = dimensions[i];
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
						dimensions.push_back(insNode);				
					
					//	Nodes currently exist at the dimension, insert as a sibling
					//		This should be inserted with nodes of the same parent
					} else {
						indTreeNode* iterateNode = dimensions[i+1];
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
	
	//Create our new node to insert
	indTreeNode* curNode = new indTreeNode;
	curNode->index = indexCounter;
	
	std::set<unsigned> tempSet = {curNode->index};
	
	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(head == nullptr){		
		head = curNode;
		indexCounter++;
		nodeCount++;
		dimensions.push_back(head);
		
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

	indTreeNode* temp = dimensions[0];
	
	
	
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
	
	
	//std::cout << "Adding neighbor..." << std::endl;
	
	//Insert into the right of the tree
	temp = dimensions[0];
	indTreeNode* ins = new indTreeNode();
	while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);
	ins->index = indexCounter;
	ins->sibling = nullptr;
	ins->parent = nullptr;
	temp->sibling = ins;
	//indexedGraph[0].push_back(std::make_pair(ins, std::make_pair(tempSet, 0)));
	
	nodeCount++;
	
	indexCounter++;
	return;

}


// Iterative function to search for a key in the tree. The function returns true
// if the key is found, else it returns false.
bool indSimplexTree::search(std::set<unsigned> simplex){
	indTreeNode* curNode = dimensions[0];
	
	for(auto i :simplex){
		while(curNode != nullptr && curNode->index != i);
		if(curNode != nullptr){ curNode = curNode->child; }
		else { return false; }
	}
	return true;
	

}


// The following function returns true if a given node has a child.
bool indSimplexTree::haveChild(indSimplexTree const* present) {

}

// A recursive function to delete a simplex from the tree.
bool indSimplexTree::deletion(indSimplexTree*& present, std::string key) {
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

double indSimplexTree::getSize(){
	//Size of node: [int + byte (*) + byte (*)] = 18 Bytes
	return nodeCount * sizeof(indTreeNode);
}

void indSimplexTree::sortAndBuildGraph(){
	std::cout << "sortAndBuildGraph" << std::endl;
	
	indTreeNode* curNode = dimensions[0];
	std::vector<graphEntry> curEntry(dimCounts[0], *(new graphEntry()));
	std::vector<std::vector<graphEntry>> ret;
	
	//Handle points first, default ordering
	while(curNode != nullptr){
		std::cout << curNode->weight << "\t" << curNode << std::endl;
		curEntry.push_back(graphEntry(curNode->simplexSet, curNode->weight, curNode));
		if(curNode->sibling == nullptr){break;}
		curNode = curNode->sibling;
	}	
	ret.push_back(curEntry);
	
	//When inserting the next few dimensions, need to first update 
	//	the indTreeNode to point to the correct graphEntry now (since we've resorted)
	//
	//	This just involves iterating the sorted vector and following the indTreeNode
	//		pointer to update the curNode; we'll use this to speed up checkFace functions

	
	
	for(int i = 1; i < dimensions.size(); i++){
		
		std::cout << "dim: " << i << std::endl;
		std::vector<graphEntry> innerEntry(dimCounts[i], *(new graphEntry()));
		curNode = dimensions[i];
		std::cout << "test" << std::endl;
		while(curNode != nullptr){
			innerEntry.push_back(graphEntry(curNode->simplexSet, curNode->weight, curNode));
			
			if(curNode->sibling == nullptr){break;}
		}		
		
		
		//SORT graph vector HERE, reassociate indTreeNode pointer
		
		
		ret.push_back(innerEntry);
	}
	
	indexedGraph = ret;
	
}


std::vector<std::vector<std::pair<std::set<unsigned>,double>>> indSimplexTree::getAllEdges(double epsilon){
	std::cout << "getAllEdges" << std::endl;
	
	if(!isSorted){
		sortAndBuildGraph();
		isSorted = true;
	}
	
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> ret;
	
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
	std::cout << "getIndexEdges" << std::endl;
	if(!isSorted){
		sortAndBuildGraph();
		isSorted = true;
	}
	return indexedGraph;
}


