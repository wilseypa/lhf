#include <string>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <iostream>
#include "simplexTree.hpp"

simplexTree::simplexTree(double _maxEpsilon, std::vector<std::vector<double>>* _distMatrix, int _maxDim){
	indexCounter = 0;
	std::cout << "simplexTree set indexCounter" << std::endl;
	distMatrix = _distMatrix;
	maxDimension = _maxDim;
	maxEpsilon = _maxEpsilon;
	simplexType = "simplexTree";
	return;
}

//**													**//
//** 				Private Functions 					**//
//**													**//
void simplexTree::recurseInsert(simplexNode* node, unsigned curIndex, int depth, double maxE, std::set<unsigned> simp){
	//Incremental insertion
	//Recurse to each child (which we'll use the parent pointer for...)
	simplexNode* temp;
	
	isSorted = false;
	double curE = 0;

	//std::cout << runningVectorIndices.size() << "\t" << runningVectorCount << "\t" << indexCounter << "\t" << distMatrix.size() << "\t" << node->index << std::endl;
	if(runningVectorIndices.size() < runningVectorCount+1){
		int offset = runningVectorCount+1 - runningVectorIndices.size();
		if(((int)node->index - offset) > distMatrix->size() || (indexCounter-offset) > distMatrix[(int)node->index - offset].size()){
			std::cout << "DistMatrix access error:" << std::endl;
			std::cout << "DistMatrix size: " << distMatrix->size() << "\tAccess Index: " << ((int)node->index - offset) << std::endl;
			std::cout << "Node Index: " << node->index << "\tOffset: " << offset << std::endl;
			std::cout << "nodeCount: " << nodeCount << "\tindexCount: " << indexCounter << std::endl;
		}
		else
			curE = (*distMatrix)[(unsigned)node->index - offset][indexCounter - offset];

	}else{
		curE = (*distMatrix)[node->index][indexCounter];
	}

	curE = curE > maxE ? curE : maxE;

	//Check if the node needs inserted at this level
	if(curE < maxEpsilon){
		simp.insert(node->index);		
		simplexNode* insNode = new simplexNode();
		insNode->index = curIndex;
		insNode->simplex = simp;
		insNode->sibling = nullptr;
		insNode->child = nullptr;
		nodeCount++;

		//Get the largest weight of this simplex
		insNode->weight = curE > node->weight ? curE : node->weight;
		maxE = insNode->weight;
		
		//if depth (i.e. 1 for first iteration) is LT weightGraphSize (starts at 1)
		if(simplexList.size() < simp.size()){
			simplexList.push_back({insNode});
		} else {
			simplexList[simp.size() - 1].insert(insNode);
		}

		//Check if the node has children already...
		if(node->child == nullptr){
			node->child = insNode;
			node->children.insert(insNode);
			insNode->parent = node;

		} else {
			temp = *node->children.rbegin();
			temp->sibling = insNode;
			insNode->parent = temp->parent;
			node->children.insert(insNode);
			temp = node->child;
			//Have to check the children now...
			if(simp.size() <= maxDimension){
				do {
					if(temp != insNode)
						recurseInsert(temp, curIndex, depth + 1, maxE, simp);
				} while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);
			}
		}
	}

	return;
}


void simplexTree::printTree(simplexNode* head){
	std::cout << "_____________________________________" << std::endl;
	if(head == nullptr){
		std::cout << "Empty tree... " << std::endl;
		return;
	}

	std::cout << "HEAD: " << head->index << "\t" << head << "\t" << head->child << "\t" << head->sibling << std::endl;
	std::cout << "sList Size: " << simplexList.size() << std::endl;


	simplexNode* current;
	
	for(int i = 0; i < simplexList.size() ; i++){
		std::cout << std::endl << "Dim: " << i << "\tSize: " << simplexList[i].size() << std::endl; 
		std::cout << "[index , address, sibling, child, parent]" << std::endl << std::endl;
		
		for(auto simplexIter = simplexList[i].begin(); simplexIter != simplexList[i].end(); simplexIter++){
			std::cout << (*simplexIter)->index << "\t" << &(*simplexIter) << "\t" << (*simplexIter)->sibling << "\t" << (*simplexIter)->child << "\t" << (*simplexIter)->parent << std::endl;
		} 

		std::cout << std::endl;
	}

	std::cout << "_____________________________________" << std::endl;
	return;
}

// Insert a node into the tree using the distance matrix and a vector index to track changes
bool simplexTree::insertIterative(std::vector<double> &currentVector, std::vector<std::vector<double>> &window){
	if(window.size() == 0){
		return true;
	}

	if(streamEval(currentVector, window)) {   // Point is deemed 'significant'
		std::vector<double> distMatrixRow = ut.nearestNeighbors(currentVector, window);

		//auto tempIndex = runningVectorIndices[0]->sibling;
		deleteIterative(runningVectorIndices[0]);
		//head = tempIndex;
		
		runningVectorIndices.erase(runningVectorIndices.begin());

		distMatrix->push_back(distMatrixRow);
		insert(distMatrixRow);

		removedSimplices++;

		return true;
	}

	return false;
}

// Insert a node into the tree using the distance matrix and a vector index to track changes
bool simplexTree::insertIterative(std::vector<double> &currentVector, std::vector<std::vector<double>> &window, int &keyToBeDeleted, int &indexToBeDeleted, std::vector<double> &distsFromCurrVec){
	if(window.size() == 0){
		return true;
	}

	if(streamEval(currentVector, window)) {   // Point is deemed 'significant'

		// deleteIterative(keyToBeDeleted);
		runningVectorIndices.erase(runningVectorIndices.begin() + indexToBeDeleted);

		insert(distsFromCurrVec);

		removedSimplices++;

		return true;
	}

	return false;
}


// Delete a node from the tree and from the distance matrix using a vector index
void simplexTree::deleteIterative(simplexNode* simplex){
	//Find what row/column of our distance matrix pertain to the vector index
	std::vector<simplexNode*>::iterator it;
	if((it = std::find(runningVectorIndices.begin(), runningVectorIndices.end(), simplex)) != runningVectorIndices.end()){

		int index = it - runningVectorIndices.begin();

		//Delete the row and column from the distance matrix based on vector index
		if(index > 0){
			
			distMatrix->erase(distMatrix->begin() + index);
		
			for(auto z : *(distMatrix)){
				if(z.size() >= index)
					z.erase(z.begin() + index);
			}
		}

		auto curNodeCount = nodeCount;

		//Delete all entries in the simplex tree with the index...
		deleteIndexRecurse(simplex->index, head);
		
		deleteWeightEdgeGraph(simplex->index);
		
		//std::cout << "Node Count reduced from " << curNodeCount << " to " << nodeCount << std::endl;
		

	} else {
		ut.writeDebug("simplexTree","Failed to find vector by index");
	}
	return;
}


void simplexTree::deleteIndexRecurse(int vectorIndex) {

    deleteIndexRecurse(vectorIndex, head);
    return;
}


void simplexTree::deleteIndexRecurse(int vectorIndex, simplexNode* curNode){
	//std::cout << "deleteIndexRecurse: " << vectorIndex << std::endl;
	if(curNode == nullptr){
		std::cout << "Empty tree" << std::endl;
		return;
	}
	
	//Handle siblings - either they need to be removed (and 'hopped') or recursed
	if(curNode->sibling != nullptr && curNode->sibling->index == vectorIndex){
		simplexNode* tempNode = curNode->sibling;
		curNode->sibling = curNode->sibling->sibling;
		deleteIndexRecurse(vectorIndex, curNode->sibling);

	} else if(curNode->sibling != nullptr){
		deleteIndexRecurse(vectorIndex, curNode->sibling);
	}

	if(curNode->index == vectorIndex){

		if(curNode == head){
			head = curNode->sibling;
			//simplexList[0].insert(simplexList.begin(), curNode->sibling); 
		}
		
		for(auto d : simplexList){
			if((*d.begin()) == curNode)
				d.insert(d.begin(), curNode->sibling);
		}
			
		deletion(curNode);
	} else if (curNode->child != nullptr && curNode->child->index == vectorIndex){
		simplexNode* tempNode = curNode->child;
		curNode->child = curNode->child->sibling;
		
		for(auto d : simplexList){
			if((*d.begin()) == tempNode)
				d.insert(d.begin(), tempNode->sibling);
		}
		
		deletion(tempNode);

	} else if(curNode->child != nullptr){
		deleteIndexRecurse(vectorIndex, curNode->child);
	}

	return;
}



// Insert a node into the tree
//
void simplexTree::insert(std::vector<double>&) {
	if(distMatrix->size() == 0){
		ut.writeDebug("simplexTree","Distance matrix is empty, skipping insertion");
		return;
	}

	//Create our new node to insert
	simplexNode* curNode = new simplexNode;
	curNode->index = indexCounter;
	std::set<unsigned> tempSet = {curNode->index};
	curNode->simplex = tempSet;
	
	//Track this index in our current window (for sliding window)
	runningVectorIndices.push_back(curNode);

	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(head == nullptr){
		root = new simplexNode;
		head = curNode;
		head->parent = root;
		root->children.insert(head);
		root->child = head;
		indexCounter++;
		runningVectorCount++;
		nodeCount++;
		simplexList.push_back({curNode});

		return;
	}

	// This needs to be a recursive span -> (or not!)
	//		if the node has a child, recurse
	//		if the node has a sibling, recurse
	//
	//root            ({!})
	//			       /
	//			      /
	//d0 -->	   | 0 | 1 | ... |
	//		        /
	//             /
	//d1 -->    | 1 | 2 | 3 |
	//           /         \
	//	        /           \
	//d2 --> | 2 | 3 |     | 4 | 5 |
	//


	//Now let's do it the new way - incremental VR;
	//	Start at the first dimension (d0);
	//		if e (d0, dn) < eps
	//			recurse down tree
	//				recurse
	//			insert to current
	//	iterate to d0->sibling

	for(auto simplexListIter = simplexList[0].begin(); simplexListIter != simplexList[0].end(); simplexListIter++){
		recurseInsert((*simplexListIter), indexCounter, 0, 0, tempSet);
	}

	//Insert into the right of the tree
	curNode->parent = root;
	root->children.insert(curNode);
	auto temp = *root->children.rbegin();
	temp->sibling = curNode;
	simplexList[0].insert(curNode);

	nodeCount++;
	indexCounter++;
	runningVectorCount++;
}

void simplexTree::deleteWeightEdgeGraph(int index){
	
	for(unsigned dim = 0; dim < simplexList.size(); dim++){
		
		for(auto simplexListIter = simplexList[dim].begin(); simplexListIter != simplexList[dim].end(); ){
			
			if((*simplexListIter)->simplex.find(index) != (*simplexListIter)->simplex.end())
				simplexList[dim].erase(*simplexListIter);
			else
				simplexListIter++;
		}
	}
	return;
}


int simplexTree::vertexCount(){
	//Return the number of vertices currently represented in the tree
	if(runningVectorIndices.size() < runningVectorCount+1)
		return runningVectorIndices.size();
	return indexCounter;
}

int simplexTree::simplexCount(){
	//Return the number of simplices currently in the tree
	return nodeCount;
}

double simplexTree::getSize(){
	//Size of node: [int + byte (*) + byte (*)] = 18 Bytes
	return nodeCount * sizeof(simplexNode);
}

//Search for a simplex from a node in the tree
simplexNode* simplexTree::find(std::set<unsigned>::iterator it, std::set<unsigned>::iterator end, simplexNode* curNode){
	simplexNode* temp = new simplexNode;
	temp->index = *it;

	while(it != end){
		auto child = curNode->children.find(temp); //Look for the sibling with the next vertex
		if(child == curNode->children.end()){ //This vertex is not in this level
			delete temp;
			return nullptr;
		} else{ //Search for the next vertex in the next level
			++it;
			temp->index = *it;
			curNode = *child;
		}
	}

	delete temp;
	return curNode;
}

std::vector<simplexNode*> simplexTree::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<simplexNode*, unsigned>& pivotPairs, bool checkEmergent){
	std::vector<simplexNode*> ret;
	simplexNode* parentNode = find(simplex.begin(), simplex.end(), root);
	if(parentNode == nullptr) return ret; //Simplex isn't in the simplex tree	

	simplexNode* tempNode;
	auto it = simplex.end();

	while(true){
		//Insert all of the children in reverse lexicographic order
		for(auto itS = parentNode->children.rbegin(); itS != parentNode->children.rend(); itS++){
			if(it == simplex.end()) ret.push_back(*itS); //All children of simplex are cofacets
			else{
				tempNode = find(it, simplex.end(), *itS); //See if cofacet is in the tree
				if(tempNode != nullptr){
					ret.push_back(tempNode);


					//If we haven't found an emergent candidate and the weight of the maximal cofacet is equal to the simplex's weight
					//		we have identified an emergent pair; at this point we can break because the interval is born and dies at the 
					//		same epsilon
					if(checkEmergent && tempNode->weight == simplexWeight){
						if(pivotPairs.find(tempNode) == pivotPairs.end()) return ret; //Check to make sure the identified cofacet isn't a pivot
						checkEmergent = false;
					}
				}
			}
		}

		//Recurse backwards up the tree and try adding vertices at each level
		--it;
		if(parentNode->parent != nullptr) parentNode = parentNode->parent;
		else break;
	}

	return ret;
}

void simplexTree::reduceComplex(){
	if(simplexList.size() == 0){
		ut.writeDebug("simplexTree","Complex is empty, skipping reduction");
		return;
	}

	//Start with the largest dimension
	ut.writeDebug("simplexTree","Reducing complex, starting simplex count: " + std::to_string(simplexCount()));
	simplexNode* cur;

	if(simplexList.size() > 0){
		for(auto i = simplexList.size()-1; i > 1; i--){

			std::vector<std::set<unsigned>> removals;
			std::vector<std::set<unsigned>> checked;

			while(checked.size() != simplexList[i].size()){
				//cur = simplexList[i];
				do {
					if(std::find(checked.begin(),checked.end(),cur->simplex) == checked.end()){
						auto ret = recurseReduce(cur, removals, checked);
						removals = ret.first;
						checked = ret.second;
					}
				} while (cur->sibling != nullptr && (cur = cur->sibling) != nullptr);
			}

			//Remove the removals
			for(auto rem : removals){
				deletion(rem);
			}

		}
	}
	ut.writeDebug("simplexTree","Finished reducing complex, reduced simplex count: " + std::to_string(simplexCount()));

	return;
}

std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> simplexTree::recurseReduce(simplexNode* simplex, std::vector<std::set<unsigned>> removals, std::vector<std::set<unsigned>> checked){
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

bool simplexTree::find(std::set<unsigned>){

	ut.writeLog("simplexTree","find(std::set<unsigned>) not implemented!");
	return 0;
}

// A recursive function to delete a simplex (and sub-branches) from the tree.
bool simplexTree::deletion(std::set<unsigned> removalEntry) {
	//Remove the entry in the simplex tree
	bool found = true;
	simplexNode* curNode = (*simplexList[0].begin());
	int d = removalEntry.size() - 1;


	/*for(int i = 0; i < dimensions[d].size(); i++){
		if(removalEntry == dimensions[d][i]->simplexSet){
			curNode = dimensions[d][i];
			dimensions[d]->erase(dimensions[d].begin() + i);
			return deletion(curNode);
		}
	}*/


	ut.writeLog("simplexTree","deletion(std::set<unsigned>) not implemented!");
	return false;

}


// A recursive function to delete a simplex (and sub-branches) from the tree.
bool simplexTree::deletion(simplexNode* removalEntry) {
	simplexNode* curNode = removalEntry;

	//Iterate to the bottom of branch in the current node
	while(curNode->child != nullptr){
		curNode->child->parent = curNode;
		curNode=curNode->child;
	}
	
	if(curNode == removalEntry && curNode->sibling != nullptr){
		deletion(curNode->sibling);
	}

	//If we did go down, remove on the way back up
	while(curNode != removalEntry){
		if(curNode->sibling != nullptr){
			deletion(curNode->sibling);
		}
		curNode = curNode->parent;

		nodeCount--;
		delete curNode->child;
		curNode->child = nullptr;
	}

	//curNode = curNode->parent;
	nodeCount--;
	delete curNode;
	return false;
}

void simplexTree::clear(){

	//Clear the simplexTree structure
	deletion(root);
	root = nullptr;

	//Clear the weighed edge graph
	for(auto z : simplexList){
		z.clear();
	}
	simplexList.clear();


	//runningVectorCount = 0;
	runningVectorIndices.clear();
	//indexCounter = 0;

}

