#include <string>
#include <vector>
#include <unistd.h>
#include <iostream>
#include "simplexTree.hpp"

simplexTree::simplexTree(double _maxEpsilon, std::vector<std::vector<double>> _distMatrix, int _maxDim){
	indexCounter = 0;
	distMatrix = _distMatrix;
	maxDim = _maxDim;
	maxEpsilon = _maxEpsilon;
	simplexType = "simplexTree";
	return;
}

void simplexTree::printTree(treeNode* head){
	
	std::cout << "_____________________________________" << std::endl;
	
	if(head == nullptr){
		std::cout << "Empty tree... " << std::endl;
		return;
	}
	
	treeNode* current;
	
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
	
	

// Insert a node into the tree
//		
void simplexTree::insert(std::vector<double>) {
	std::cout << "Insert Node - Index: " << indexCounter << "\tNode Count: " << nodeCount << "\tSize: " << getSize() << std::endl;
	
	//printTree(head);
	
	
	//Create our new node to insert
	treeNode* curNode = new treeNode;
	curNode->index = indexCounter;
	
	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(head == nullptr){		
		head = curNode;
		indexCounter++;
		nodeCount++;
		std::cout << "HEAD: " << head << std::endl;
		dimensions.push_back(head);
		
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

	
	//Loop each dimensional list; start at dn, then dn-1, ..., d0
	for(int i = (dimensions.size()-1 < maxDim ? dimensions.size() - 1 : maxDim); i >= 0; i--){
		
		
		
		//Loop each sibling in the current dimensional list (as a do-while, break when no more siblings)
		curNode = dimensions[i];
		do{
			//Determine if the index needs to inserted below the current node
			//	First, check the distance matrix for curNodeIndex v. curIndex
			//	Second (if first is true), check all parents dist matrix of parentNodeIndex v. curIndex
			//astd::cout << "TEST: " << i << std::endl;
			if(distMatrix[curNode->index][indexCounter] < maxEpsilon){
				//	This node is a candidate, now we need to check each parent of the current branch
				//		to ensure the distance matrix entry for parentNodeindex v. curIndex is < maxEpsilon

				bool ins = true;
				treeNode* parentNode = curNode; //Temporarily store the parent so we don't lose track of our current
			
				// This loop will check each parent and set whether to insert a node
				do{
					if(distMatrix[parentNode->index][indexCounter] > maxEpsilon){
						ins = false;
					}	
				}while(ins && parentNode->parent != nullptr && (parentNode = parentNode->parent) != nullptr);
				
				//Insert the node; the distance to each parent is less than epsilon
				if(ins){
					//Allocate a new node to be inserted into the tree
					treeNode* insNode = new treeNode;
					
					//This new node's parent is the current node we're indexing on
					insNode->parent = curNode;
					insNode->index = indexCounter;
					insNode-> sibling = nullptr;
				
					//Check the dimensional list
					//	if no nodes exist at the dimension, this is the first
					if(i == dimensions.size() - 1){
						dimensions.push_back(insNode);				
					
					//	Nodes currently exist at the dimension, insert as a sibling
					//		This should be inserted with nodes of the same parent
					} else {
						treeNode* iterateNode = dimensions[i+1];
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
	
	
	
	
	//std::cout << "Adding neighbor..." << std::endl;
	
	//Insert into the right of the tree
	treeNode* temp = dimensions[0];
	treeNode* ins = new treeNode();
	while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);
	ins->index = indexCounter;
	ins->sibling = nullptr;
	ins->parent = nullptr;
	temp->sibling = ins;
	
	
	nodeCount++;
	
	indexCounter++;
	return;

}


// Iterative function to search for a key in the tree. The function returns true
// if the key is found, else it returns false.
bool simplexTree::search(std::string key){


}


// The following function returns true if a given node has a child.
bool simplexTree::haveChild(simplexTree const* present) {

}

// A recursive function to delete a simplex from the tree.
bool simplexTree::deletion(simplexTree*& present, std::string key) {
	return false;
}

int simplexTree::vertexCount(){
	//Return the number of vertices in the tree
	return nodeCount;
}

int simplexTree::simplexCount(){
	//Return the number of simplices in the tree
	
	return -1;
}

double simplexTree::getSize(){
	//Size of node: [double + double + byte (*) + byte (*)] = 18 Bytes
	return nodeCount * 18;
}
