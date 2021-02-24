#include <string>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <iostream>
// #include <typeinfo>
#include "simplexTree.hpp"

simplexTree::simplexTree(double _maxEpsilon, int _maxDim){
	indexCounter = 0;
	maxDimension = _maxDim;
	maxEpsilon = _maxEpsilon;
	simplexType = "simplexTree";
	return;
}

void simplexTree::outputComplex(){
	return printTree(root->child);
}

void simplexTree::recurseInsert(simplexTreeNode* node, unsigned curIndex, int depth, double maxE, std::set<unsigned> simp){
	//Incremental insertion
	//Recurse to each child (which we'll use the parent pointer for...)
	simplexTreeNode* temp;
	double curE = 0;


	if(runningVectorIndices.size() < runningVectorCount){

		//Get the positions of the vector in the runningVectorIndices array
		auto nodeIndex = std::find(runningVectorIndices.begin(), runningVectorIndices.end(), node->simpNode->index);

		if((std::distance(runningVectorIndices.begin(), nodeIndex)) > distMatrix->size() || (indexCounter - (runningVectorCount - 1)) > (*distMatrix)[std::distance(runningVectorIndices.begin(), nodeIndex)].size()){
			std::cout << "DistMatrix access error:" << std::endl;
			std::cout << "\tAttempting to access distMatrix indexes: " << node->simpNode->index << " x " << indexCounter << std::endl;
			std::cout << "\tDistMatrix size: " << (*distMatrix).size() << std::endl;
			std::cout << "\trviCount: " << runningVectorCount << "\t rviSize: " << runningVectorIndices.size() << "\tOffset: " << simplexOffset << "\tIC: " << indexCounter << std::endl;
			std::cout << "\tOffset Indices: " << node->simpNode->index - (runningVectorCount - 1) << " x " << indexCounter - (runningVectorCount - 1) << std::endl;
			std::cout << "\tBackwards size: " << distMatrix[indexCounter - (runningVectorCount - 1)].size() << std::endl;
			std::cout << "\tRow Size: " << distMatrix[indexCounter - (runningVectorCount - 1)].size()  << "\tCurIndex: " << curIndex << std::endl;
			std::cout << "\tNode Index: " << std::distance(runningVectorIndices.begin(), nodeIndex) << std::endl;
		}
		else {
		    curE = *((*distMatrix)[std::distance(runningVectorIndices.begin(), nodeIndex)].rbegin());
		}


	} else {
		curE = (*distMatrix)[node->simpNode->index][indexCounter];
	}


	curE = curE > maxE ? curE : maxE;

	//Check if the node needs inserted at this level
	if(curE <= maxEpsilon){
		simp.insert(node->simpNode->index);
		//Get the largest weight of this simplex
		maxE = curE > node->simpNode->weight ? curE : node->simpNode->weight;
		simplexTreeNode* insNode = new simplexTreeNode(simp, maxE);
		insNode->simpNode->index = curIndex;
	    insNode->simpNode->hash = nodeCount;
		nodeCount++;


		//if depth (i.e. 1 for first iteration) is LT weightGraphSize (starts at 1)
		//if(simplexList.size() < simp.size())
		//	simplexList.push_back({});
		//simplexList[simp.size() - 1].insert(insNode->simp);

		//Check if the node has children already...
		if(node->child == nullptr){
			node->child = insNode;
			insNode->parent = node;
			//node->children.insert(insNode);

		//Node has children, add to the beginning of children
		} else {
			//Add new node as the first child. Shift existing children to left.
			insNode->parent = node;
			insNode->sibling = node->child;
			node->child = insNode;
			//node->children.insert(insNode);

			temp = insNode->sibling;

			//Have to check the children now...
			if(simp.size() <= maxDimension){
				do {
					recurseInsert(temp, curIndex, depth + 1, maxE, simp);
				} while((temp = temp->sibling) != nullptr);
			}
		}
	}

	return;
}


void simplexTree::printTree(simplexTreeNode* headPointer){
	std::cout << "_____________________________________" << std::endl;
	if(root->child == nullptr){
		std::cout << "Empty tree... " << std::endl;
		return;
	}
	auto temp = headPointer;

	std::cout << "ROOT: " << headPointer->simpNode->index << "\t" << headPointer << "\t" << headPointer->child << "\t" << headPointer->sibling << std::endl;
	

	std::cout << "[index , address, sibling, child, parent]" << std::endl << std::endl;

	for(auto simplexIter = headPointer; simplexIter != nullptr; simplexIter = simplexIter->sibling){
		std::cout << simplexIter->simpNode->index << "\t";
		std::cout << simplexIter << "\t" ;
		std::cout << simplexIter->sibling << "\t" ;
		std::cout << simplexIter->child << "\t" ;
		std::cout << simplexIter->parent << "\t";
		ut.print1DVector(simplexIter->simpNode->simplex);
		
		//temp = simplexIter;
	}

	std::cout << "_____________________________________" << std::endl;
	
	std::cout << "Children of root->child (" << temp->sibling->sibling->sibling->sibling->child << ")" << std::endl << std::endl;
	
	for(auto simplexIter = temp->sibling->sibling->sibling->sibling->child; simplexIter != nullptr; simplexIter = simplexIter->sibling){
		std::cout << simplexIter->simpNode->index << "\t";
		std::cout << simplexIter << "\t" ;
		std::cout << simplexIter->sibling << "\t" ;
		std::cout << simplexIter->child << "\t" ;
		std::cout << simplexIter->parent << "\t";
		ut.print1DVector(simplexIter->simpNode->simplex);
	}
	
	return;
}

// Insert a node into the tree using the distance matrix and a vector index to track changes
bool simplexTree::insertIterative(std::vector<double> &currentVector, std::vector<std::vector<double>> &window){
	if(window.size() == 0){
		return true;
	}

	if(streamEval(currentVector, window)) {   // Point is deemed 'significant'

		//Delete the oldest point in the window
		deleteIterative(runningVectorIndices[0]);
		runningVectorIndices.erase(runningVectorIndices.begin());

		//Create distance matrix row of current vector to each point in the window
		std::vector<double> distsCurrVec = ut.nearestNeighbors(currentVector, window);
		// std::vector<double> distMatrixRow = ut.nearestNeighbors(currentVector, window);

		distsCurrVec.erase(distsCurrVec.begin());

		//Insert the new point into the distance matrix and complex
		for(int i = 0; i < (*distMatrix).size(); i++) {
            (*distMatrix)[i].push_back(distsCurrVec[i]);
		}

		distsCurrVec.push_back(0);
//
//		std::vector<double> distMatLastRow( window.size() );
		distMatrix->push_back(distsCurrVec);

//        distMatrix->push_back(distMatrixRow);

		insert();

		removedSimplices++;

		return true;
	}

	return false;
}

// Insert a node into the tree using the distance matrix and a vector index to track changes
bool simplexTree::insertIterative(std::vector<double> &currentVector, std::vector<std::vector<double>> &window, int &keyToBeDeleted, int &indexToBeDeleted){
	if(window.size() == 0){
		return true;
	}

	if(streamEval(currentVector, window)) {   // Point is deemed 'significant'

//	     std::cout << "========================== Simplex tree after insertion ==========================" << '\n';
//		 printTree(root);

		std::cout << "indexToBeDeleted = " << indexToBeDeleted << '\n';
		deleteIndexRecurse(keyToBeDeleted);
		runningVectorIndices.erase(runningVectorIndices.begin() + indexToBeDeleted);

//		 std::cout << "========================== Simplex tree after deletion ==========================" << '\n';
//		 printTree(root);

		insert();

		removedSimplices++;

		return true;
	}

	return false;
}


// Delete a node from the tree and from the distance matrix using a vector index
void simplexTree::deleteIterative(int simplexIndex){

	//Find what row/column of our distance matrix pertain to the vector index
	std::vector<int>::iterator it;
	if((it = std::find(runningVectorIndices.begin(), runningVectorIndices.end(), simplexIndex)) != runningVectorIndices.end()){

		//Index holds the index in the runningVectorIndices array of the simplexNode pointer
		int index = std::distance(runningVectorIndices.begin(), it);
		std::cout << "index = " << index << '\n';

		//Delete the row and column from the distance matrix based on vector index
		//	This corresponds to the index into the runnningVectorIndices array

		//Delete Row[index]
		distMatrix->erase(distMatrix->begin() + index);

		//Delete column[index] (row[][index])
		for(int i = 0; i < (*distMatrix).size(); i++) {
            if((*distMatrix)[i].size() >= index)
                (*distMatrix)[i].erase((*distMatrix)[i].begin() + index);
		}


        auto curNodeCount = nodeCount;

		//Delete all entries in the simplex tree
		//printTree(root);
		// std::cout << "simplex->index " << simplex->index << '\n';
		deleteIndexRecurse(simplexIndex);

		//printTree(root);

	} else {
		ut.writeDebug("simplexTree","Failed to find vector by index");
	}
	return;
}


void simplexTree::deleteIndexRecurse(int vectorIndex) {
    std::cout << "deleteIndexRecurse vectorIndex = " << vectorIndex << '\n';

    // Since child index is always higher than parent index, no "top node" between root->child and vectorIndex
    // can contain a subtree that has a node->index =  vectorIndex. Therefore, skip to the top node whose sibling
    // has vectorIndex.
    simplexTreeNode* curNode = root->child;

    if(curNode->sibling != nullptr)
    {
        while(curNode->sibling->simpNode->index > vectorIndex)
        {
            curNode = curNode->sibling;
        }
    }

    deleteIndexRecurse(vectorIndex, curNode);
    return;
}


void simplexTree::deleteIndexRecurse(int vectorIndex, simplexTreeNode* curNode){
	if(curNode == nullptr){
		std::cout << "Empty tree" << std::endl;
		return;
	}

	//Handle siblings - either they need to be removed (and 'hopped') or recursed
	//	Only need to recurse if the vector could be a member of tree (i.e. < vectorIndex)
	if(curNode->sibling != nullptr && curNode->sibling->simpNode->index == vectorIndex){   // Current node's sibling is to be deleted

		//Map the current node's sibling to the node following the node for deletion
		simplexTreeNode* tempNode = curNode->sibling;
		curNode->sibling = curNode->sibling->sibling;

		//Delete the orphaned node now that sibling remapping is handled
		//	The sibling still points to the following node, so this is valid
		deleteIndexRecurse(vectorIndex, tempNode);

	} else if(curNode->sibling != nullptr){
		//Recurse to the next sibling and check for deletion
		deleteIndexRecurse(vectorIndex, curNode->sibling);
	}


	//Handle self; if this is the vector index delete branches/nodes
	//	Assume all children / sibling pointers have been alleviated in calling function
	if(curNode->simpNode->index == vectorIndex){
		//Ensure we aren't the first node of the tree
		if(curNode == root->child){
			root->child = curNode->sibling;
		}

		//Remove our index from the parent
		//curNode->parent->children.erase(curNode);

		//Delete this node and sub-nodes in the tree
		deletion(curNode);

	//If this isn't the vectorIndex, need to look at children and remove or recurse
	//	Only need to recurse if the vector could be a member of tree (i.e. < vectorIndex)
	} else if (curNode->child != nullptr && curNode->child->simpNode->index == vectorIndex){
		simplexTreeNode* tempNode = curNode->child;
		curNode->child = curNode->child->sibling;

//		for(auto d : simplexList){
//			if((*d.begin()) == tempNode) {
//                d.insert(d.begin(), tempNode->sibling);
//			}
//		}

		deleteIndexRecurse(vectorIndex, tempNode);

	} else if( curNode->child != nullptr && curNode->child->simpNode->index > vectorIndex ) {
	    // If curNode->child->index < vectorIndex, the subtree rooted at curNode->child cannot contain vectorIndex.
		deleteIndexRecurse(vectorIndex, curNode->child);

	}

	return;
}


// Insert a node into the tree
//
void simplexTree::insert() {
	if(distMatrix->size() == 0){
		ut.writeDebug("simplexTree","Distance matrix is empty, skipping insertion");
		return;
	}
	
	//Create our new node to insert (Ref Count = 1)
	simplexTreeNode* insNode = new simplexTreeNode({(unsigned)indexCounter}, 0); 	
	insNode->simpNode->index = indexCounter;

	//Track this index in our current window (for sliding window)
	runningVectorIndices.push_back(insNode->simpNode->index);

	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(root == nullptr){
		root = new simplexTreeNode();
		insNode->parent = root;
		root->child = insNode;
		//root->children.insert(insNode);
		indexCounter++;
		runningVectorCount++;
		nodeCount++;

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

	runningVectorCount++;

	for(auto it = root->child; it != nullptr; it = it->sibling){
		recurseInsert(it, indexCounter, 0, 0, {(unsigned)indexCounter});
	}

	//Insert into the right of the tree
	//Child points to last child, and sibling points backwards
	insNode->parent = root;
	insNode->sibling = root->child;
	root->child = insNode;
	//root->children.insert(insNode);

	//simplexList[0].insert(insNode);

	// runningVectorCount++;
	insNode->simpNode->hash = nodeCount;
	nodeCount++;
	indexCounter++;
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
simplexTree::simplexTreeNode* simplexTree::find(std::set<unsigned>::iterator begin, std::set<unsigned>::iterator end, simplexTreeNode* curNode){
	auto it = begin;

	while(it != end){
		simplexTreeNode* ptr;
		for(ptr = curNode->child; ptr != nullptr; ptr = ptr->sibling){
			if(ptr->simpNode->index == (*it)){
				++it;
				curNode = ptr;
				break;
			}
		}
		if(ptr == nullptr) return nullptr;
	}
	
	return curNode;
}


std::vector<std::set<simplexNode_P, cmpByWeight>> simplexTree::getAllEdges(){
	std::vector<std::set<simplexNode_P, cmpByWeight>> ret(maxDimension + 1, std::set<simplexNode_P, cmpByWeight>());
	if(root != nullptr)
		recurseGetEdges(ret, root, 0, maxDimension);
	return ret;
}

void simplexTree::recurseGetEdges(std::vector<std::set<simplexNode_P, cmpByWeight>> &edgeList, simplexTreeNode* current, int depth, int maxDepth){
	for(auto ptr = current->child; ptr != nullptr; ptr = ptr->sibling){
		
		edgeList[depth].insert(ptr->simpNode);
		
		if(ptr->child != nullptr && depth+1 <= maxDepth)
			recurseGetEdges(edgeList, ptr, depth+1, maxDepth);
	
	}
	return;
}

std::vector<simplexNode*> simplexTree::getAllCofacets(simplexNode_P simp){
	return getAllCofacets(simp, std::unordered_map<long long, simplexNode_P>(), false);
}

std::vector<simplexNode*> simplexTree::getAllCofacets(simplexNode_P simp, const std::unordered_map<long long, simplexNode_P>& pivotPairs, bool checkEmergent){
	std::vector<simplexNode*> ret;
	simplexTreeNode* parentNode = find(simp->simplex.begin(), simp->simplex.end(), root);
	if(parentNode == nullptr) return ret; //Simplex isn't in the simplex tree

	simplexTreeNode* tempNode;
	auto it = simp->simplex.end();

	while(true){
		//Insert all of the children in reverse lexicographic order
		for(auto ptr = parentNode->child; ptr != nullptr; ptr = ptr->sibling){
			if(it == simp->simplex.end()) {
				ret.push_back(ptr->simpNode.get()); //All children of simplex are cofacets
			} else {
				tempNode = find(it, simp->simplex.end(), ptr); //See if cofacet is in the tree
				if(tempNode != nullptr){
					ret.push_back(tempNode->simpNode.get());
					
					//If we haven't found an emergent candidate and the weight of the maximal cofacet is equal to the simplex's weight
					//		we have identified an emergent pair; at this point we can break because the interval is born and dies at the
					//		same epsilon
					if(checkEmergent && tempNode->simpNode->weight == simp->weight){
						if(pivotPairs.find(tempNode->simpNode->hash) == pivotPairs.end()) return ret; //Check to make sure the identified cofacet isn't a pivot
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

std::vector<simplexNode_P> simplexTree::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<simplexNode_P, simplexNode_P>& pivotPairs, bool checkEmergent){
	std::vector<simplexNode_P> ret;
	simplexTreeNode* parentNode = find(simplex.begin(), simplex.end(), root);
	if(parentNode == nullptr) return ret; //Simplex isn't in the simplex tree

	simplexTreeNode* tempNode;
	auto it = simplex.end();

	while(true){
		//Insert all of the children in reverse lexicographic order
		for(auto ptr = parentNode->child; ptr != nullptr; ptr = ptr->sibling){
			if(it == simplex.end()) {
				ret.push_back(ptr->simpNode); //All children of simplex are cofacets
			} else {
				tempNode = find(it, simplex.end(), ptr); //See if cofacet is in the tree
				if(tempNode != nullptr){
					ret.push_back(tempNode->simpNode);
					
					//If we haven't found an emergent candidate and the weight of the maximal cofacet is equal to the simplex's weight
					//		we have identified an emergent pair; at this point we can break because the interval is born and dies at the
					//		same epsilon
					if(checkEmergent && tempNode->simpNode->weight == simplexWeight){
						if(pivotPairs.find(tempNode->simpNode) == pivotPairs.end()) return ret; //Check to make sure the identified cofacet isn't a pivot
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
	simplexTreeNode* cur;

	if(simplexList.size() > 0){
		for(auto i = simplexList.size()-1; i > 1; i--){

			std::vector<std::set<unsigned>> removals;
			std::vector<std::set<unsigned>> checked;

			while(checked.size() != simplexList[i].size()){
				//cur = simplexList[i];
				do {
					if(std::find(checked.begin(),checked.end(),cur->simpNode->simplex) == checked.end()){
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

std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> simplexTree::recurseReduce(simplexTreeNode* simplex, std::vector<std::set<unsigned>> removals, std::vector<std::set<unsigned>> checked){
	checked.push_back(simplex->simpNode->simplex);
	auto subsets = ut.getSubsets(simplex->simpNode->simplex);
	std::set<unsigned> maxFace;

	bool canRemove = true;

	//Look at each face
	/*for(auto face : subsets){

		//Check if the face is shared; if so, recurse
		for(auto simp : simplexList[simplex->simpNode->simplex.size() - 1]){

			if(simp.get() != simplex && std::find(checked.begin(), checked.end(), simp->simpNode->simplex) == checked.end()){
				auto sDiff = ut.symmetricDiff(simp->simpNode->simplex, face,true);

				//This point intersects;
				if (sDiff.size() == 1){
					auto ret = recurseReduce(simp.get(), removals, checked);
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
	}*/

	return std::make_pair(removals, checked);

}

bool simplexTree::find(std::set<unsigned>){

	ut.writeLog("simplexTree","find(std::set<unsigned>) not implemented!");
	return 0;
}

// A recursive function to delete a simplex (and sub-branches) from the tree.
bool simplexTree::deletion(std::set<unsigned> removalEntry) {
	//Remove the entry in the simplex tree
	// bool found = true;
	// simplexNode_P curNode = (*simplexList[0].begin());
	// int d = removalEntry.size() - 1;


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
bool simplexTree::deletion(simplexTreeNode* removalEntry) {
	simplexTreeNode* curNode = removalEntry;

	//Iterate each child and delete
	/*auto it = curNode->children.begin();
	while(it != curNode->children.end()){
		deletion(*(it));
		it++;
	}
	nodeCount--;*/

	//Remove from the simplex list
	//simplexList[curNode->simplex.size() - 1].erase(simplexList[curNode->simplex.size() - 1].find(curNode));

	//delete curNode;
	return false;
}

void simplexTree::clear(){
	//Clear the simplexTree structure
	head = nullptr;
	root = nullptr;

	simplexList.clear();

	simplexOffset = runningVectorCount;
	runningVectorIndices.clear();
	runningVectorCount = 0;
	indexCounter = 0;
	nodeCount = 0;

	return;

}

