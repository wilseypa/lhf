#include <string>
#include <vector>
#include <iostream>
#include "simplexTree.hpp"

simplexTree::simplexTree(std::vector<std::vector<double>> _distMatrix){
	indexCounter = 0;
	distMatrix = _distMatrix;
	return;
}


// simplexTree constructor
simplexTree::simplexTree(double _maxEpsilon, std::vector<std::vector<double>> _distMatrix){
	indexCounter = 0;
	distMatrix = _distMatrix;
	
	simplexType="simplexTree";
	maxEpsilon = _maxEpsilon;
	
}

void simplexTree::recurse(treeNode* node, int curIndex){
	//If this node has a child, recurse children first
	
	//std::cout << "Check child : " << node->index << std::endl;
	/*if( node->child != nullptr ){
		recurse(node->child, curIndex);
	}

	//std::cout << "Check neighbors: " << node->index << std::endl;
	//If this node has neighbors (next), recurse neighbors
	if( node->sibling != nullptr){
		recurse(node->sibling, curIndex);
	}

	//Check if this node needs insertion, if so, make a copy and insert
	if(distMatrix[node->index][curIndex] < maxEpsilon){		
		
		treeNode* insNode = new treeNode;
		treeNode* temp;
		
		//std::cout << "\tCreating... " << std::endl;
		
		insNode->index = curIndex;
		
		//std::cout << "\tCreated new child node..." << std::endl;
		
		if(node->child == nullptr && checkParent(node, curIndex)){
			node->child = insNode;
			insNode->parent = node;
			nodeCount++;
		} else if (checkParent(node,curIndex)) {
			// Children already exist; iterate until we find empty slot
			//temp = node->child;
			while(temp->sibling != nullptr){temp = temp->sibling;};
			temp->sibling = insNode;	
			insNode->parent = temp->parent;	
			nodeCount++;			
		}	
	}*/
	
	return;
}


bool simplexTree::checkParent(treeNode* node, int curIndex){
	bool retVal = true;
	
	if(node->parent != nullptr)
		retVal = checkParent(node->parent, curIndex);
	
	if(distMatrix[node->index][curIndex] >= maxEpsilon)
		retVal = false;
		
	return retVal;
	
}

// Insert a node into the tree
//		
void simplexTree::insert(std::vector<double>) {
	std::cout << "Insert Node - Index: " << indexCounter << "\tNode Count: " << nodeCount << "\tSize: " << getSize() << std::endl;
	
	
	//Create our new node to insert
	treeNode* curNode = new treeNode;
	curNode->index = indexCounter;
	
	//std::cout << "Created Node" << std::endl;
	
	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(head == nullptr){		
		head = curNode;
		indexCounter++;
		nodeCount++;
		std::cout << "HEAD: " << head << std::endl;
		dimensions.push_back(head);
		
		return;
	} else if (dimensions.size() ==1){
		head->child = curNode;
		curNode->parent = head;
		
		std::cout << "HEAD: " << head << std::endl;
		std::cout << "CurNode: " << curNode << std::endl;
		std::cout << "curNodeParent: " << curNode->parent << std::endl;
		indexCounter++;
		nodeCount++;
		dimensions.push_back(curNode);
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
	

	// For each node, check if the current node is inserted anywhere
		
	//std::cout << "recursing... " << head->index << std::endl;	
	//recurse(head, indexCounter);
	
	std::cout << "Test: " << dimensions.size() << std::endl;
	
	//Loop each dimensional list
	for(int i = dimensions.size() - 1; i > 0; i--){
		std::cout << "a " << i << "\t" << dimensions[i] << std::endl;
		//Loop each sibling in the dimensional list
		
		treeNode* prevNode = nullptr;
		curNode = dimensions[i];
		do{
			
			std::cout << "b" << std::endl;
			//Determine if the index needs to inserted below the current node
			bool ins = true;
			treeNode* parentNode = curNode;
			
			do{
				if(distMatrix[parentNode->index][indexCounter] > maxEpsilon){
					ins = false;
					std::cout << "weight > maxE" << std::endl;
					break;
				}	
				parentNode = parentNode->parent;
				
			}while(ins && parentNode != head);
			std::cout << "c" << std::endl;
			
			if(ins){
				//Add the node to the list
				treeNode* insNode = new treeNode;
				insNode->parent = curNode;
				insNode->index = indexCounter;
				
				std::cout << "d" << std::endl;
				//If no nodes currently exist at the list level
				if(i == dimensions.size() - 1){
					dimensions.push_back(insNode);
					curNode->child = insNode;				
				} else if(curNode->child == nullptr){
					std::cout << "e" << std::endl;
					//Current node has no children
					curNode->child = insNode;
					
					std::cout << "f" << std::endl;
					//Insert into the previous subnode
					if(prevNode != nullptr){
						
						prevNode = prevNode->child;
						std::cout << "g2" << std::endl;
						while(prevNode->sibling != nullptr && (prevNode = prevNode->sibling) != nullptr){
							//std::cout << prevNode << "\t";
						}
					}
					std::cout << "g1" << std::endl;
					prevNode->sibling = insNode;
					
					//Add the next subnode
					if((curNode = curNode->sibling) != nullptr){
						insNode->sibling = curNode->child;				
					}
				} else {
					treeNode* temp = curNode->child;
					treeNode* temp_prev = temp;
					while( temp->sibling != nullptr && (temp_prev = temp)  &&  (temp = temp->sibling) && temp->sibling->parent == curNode);
					std::cout << "e2" << std::endl;
					temp_prev->sibling = insNode;
					insNode->sibling = temp;
				}
			}			
			
			prevNode = curNode;
		}while(curNode->sibling != nullptr && (curNode = curNode->sibling));
	}
	
	
	
	
	//std::cout << "Adding neighbor..." << std::endl;
	
	//Insert into the right of the tree
	treeNode* temp = head;
	while(temp->sibling != nullptr){temp = temp->sibling;};
	temp->sibling = curNode;	
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
