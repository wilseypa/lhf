#include <string>
#include <vector>
#include <iostream>
#include "simplexTree.hpp"

simplexTree::simplexTree(){return;}

// simplexTree constructor
simplexTree::simplexTree(double maxEpsilon){
	simplexType="simplexTree";
	maxEpsilon = maxEpsilon;
	this->isLeaf = false;

	for (int i = 0; i < MAX_POINTS; i++)
		this->character[i] = nullptr;
}

double simplexTree::getSize(){
	size_t size = 0;
	
	//Calculate size of original data
	
	return size;
}


// Iterative function to insert a simplex, i.e. a key, in the tree.
void simplexTree::insert(std::string key) {

	// Start from the root.
	simplexTree* present = this;
	for (int i = 0; i < key.length(); i++) {
		// Create a new node if it doesn't exist.
		if (present->character[key[i]] == nullptr)
			present->character[key[i]] = new simplexTree(maxEpsilon);

		// Go to the next node.
		present = present->character[key[i]];
	}
	
	// Set the current node as a leaf.
	present->isLeaf = true;

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
	
	return -1;
}

int simplexTree::simplexCount(){
	//Return the number of simplices in the tree
	
	return -1;
}
