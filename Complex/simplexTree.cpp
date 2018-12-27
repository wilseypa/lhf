#include <string>
#include <vector>
#include <iostream>
#include "simplexTree.hpp"

simplexTree::simplexTree(std::vector<std::vector<double>> _distMatrix){
	head = tail;
	indexCounter = 0;
	distMatrix = _distMatrix;
	return;
}


// simplexTree constructor
simplexTree::simplexTree(double _maxEpsilon, std::vector<std::vector<double>> _distMatrix){
	head = tail;
	indexCounter = 0;
	distMatrix = _distMatrix;
	
	simplexType="simplexTree";
	maxEpsilon = _maxEpsilon;
	
}

void simplexTree::recurse(treeNode* node, int curIndex){
	//If this node has a child, recurse children first
	
	//std::cout << "Check child : " << node->index << std::endl;
	if( node->child != nullptr ){
		recurse(node->child, curIndex);
	}

	//std::cout << "Check neighbors: " << node->index << std::endl;
	//If this node has neighbors (next), recurse neighbors
	if( node->next != nullptr ){
		recurse(node->next, curIndex);
	}

	//std::cout << "Inserting new node... : " << maxEpsilon << " [ " << node->index << " , " << curIndex << " ] "<< std::endl;
	//std::cout << "distMatrix: " << distMatrix[node->index][curIndex] << std::endl;
	
	//Check if this node needs insertion, if so, make a copy and insert
	if(distMatrix[node->index][curIndex] < maxEpsilon){
		
		treeNode* insNode = new treeNode;
		treeNode* temp;
		
		//std::cout << "\tCreating... " << std::endl;
		
		insNode->index = curIndex;
		insNode->weight = distMatrix[node->index][curIndex];
		
		//std::cout << "\tCreated new child node..." << std::endl;
		
		if(node->child == nullptr){
			node->child = insNode;
			nodeCount++;
		} else {
			// Children already exist; iterate until we find empty slot
			temp = node->child;
			while(temp->next != nullptr){temp = temp->next;};
			temp->next = insNode;			
			nodeCount++;			
		}	
	}
	
	return;
}

// Insert a node into the tree
//		
void simplexTree::insert(std::vector<double>) {
	//std::cout << "Inserting node: " << indexCounter << std::endl;
	
	//Create our new node to insert
	treeNode* curNode = new treeNode;
	curNode->index = indexCounter;
	curNode->weight = 0;
	
	//std::cout << "Created Node" << std::endl;
	
	//Check if this is the first node (i.e. head==tail)
	//	If so, initialize the head and tail nodes
	if(head == nullptr){
		//std::cout << "First node!" << std::endl;
		head = curNode;
		indexCounter++;
		nodeCount++;
		return;
	}
		
	// This needs to be a recursive span ->
	//		if the node has a child, recurse
	//		if the node has a sibling, recurse
	//
	//			| 0 | 1 | ... |
	//			/
	//         /
	//		| 1 | 2 | 3 | 
	//      /          \
	//	   /            \
	//   | 2 | 3 |     | 4 | 5 |

	// For each node, check if the current node is inserted anywhere
		
	//std::cout << "recursing... " << head->index << std::endl;	
		
	recurse(head, indexCounter);
	
	//std::cout << "Adding neighbor..." << std::endl;
	
	//Insert into the right of the tree
	treeNode* temp = head;
	while(temp->next != nullptr){temp = temp->next;};
	temp->next = curNode;	
	nodeCount++;
	
	indexCounter++;
	return;

}


// Iterative function to search for a key in the tree. The function returns true
// if the key is found, else it returns false.
bool simplexTree::search(std::string key){
	
	// Return false if the tree is empty.
	if (this == nullptr)
		return false;

	simplexTree* present = this;

	for (int i = 0; i < key.length(); i++) {
		// Go to the next node.
		present = present->character[key[i]];

		// If the key queried for is invalid, i.e. if the end of the tree is reached
		if (present == nullptr)
			return false;
	}

	// If the current node is a leaf and we have reached the end of the string, return true.
	return present->isLeaf;

}


// The following function returns true if a given node has a child.
bool simplexTree::haveChild(simplexTree const* present) {
	for (int i = 0; i < MAX_POINTS; i++)
		if (present->character[i])
			return true;    // Child found.

	return false;
}

// A recursive function to delete a simplex from the tree.
bool simplexTree::deletion(simplexTree*& present, std::string key) {

	// Return if the tree is empty
	if (present == nullptr)
		return false;

	// If the end of the key is not reached.
	if (key.length()) {
		// recursively search for the node corresponding to the next character
		// in the key and if it returns true, delete current node given it is a non-leaf node.
		if (present != nullptr &&
				present->character[key[0]] != nullptr &&
				deletion(present->character[key[0]], key.substr(1)) &&
				present->isLeaf == false) {

			if (!haveChild(present)) {
				delete present;
				present = nullptr;
				return true;
			} else {
				return false;
			}
		}
	}

	// If the end of the key is reached
	if (key.length() == 0 && present->isLeaf) {
		// If the current node is a leaf node and doesn't have any child
		if (!haveChild(present)) {
			// Delete the current node.
			delete present;
			present = nullptr;

			// Delete the non-leaf parent nodes.
			return true;
		} else {
			// Mark current node as a non-leaf node (it is not deleted).
			present->isLeaf = false;

			// Its parent nodes are not deleted.
			return false;
		}
	}

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
