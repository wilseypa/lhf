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

void simplexTree::outputComplex(){
	return printTree(root->child);
}

void simplexTree::recurseInsert(simplexNode* node, unsigned curIndex, int depth, double maxE, std::set<unsigned> simp){
	//Incremental insertion
	//Recurse to each child (which we'll use the parent pointer for...)
	simplexNode* temp;
	double curE = 0;

	if(runningVectorIndices.size() < runningVectorCount){

		//Get the positions of the vector in the runningVectorIndices array
		auto nodeIndex = std::find(runningVectorIndices.begin(), runningVectorIndices.end(), node);

		if((std::distance(runningVectorIndices.begin(), nodeIndex)) > distMatrix->size() || (indexCounter - (runningVectorCount - 1)) > (*distMatrix)[std::distance(runningVectorIndices.begin(), nodeIndex)].size()){
			std::cout << "DistMatrix access error:" << std::endl;
			std::cout << "\tAttempting to access distMatrix indexes: " << node->index << " x " << indexCounter << std::endl;
			std::cout << "\tDistMatrix size: " << (*distMatrix).size() << std::endl;
			std::cout << "\trviCount: " << runningVectorCount << "\t rviSize: " << runningVectorIndices.size() << "\tOffset: " << simplexOffset << "\tIC: " << indexCounter << std::endl;
			std::cout << "\tOffset Indices: " << node->index - (runningVectorCount - 1) << " x " << indexCounter - (runningVectorCount - 1) << std::endl;
			std::cout << "\tBackwards size: " << distMatrix[indexCounter - (runningVectorCount - 1)].size() << std::endl;
			std::cout << "\tRow Size: " << distMatrix[indexCounter - (runningVectorCount - 1)].size()  << "\tCurIndex: " << curIndex << std::endl;
			std::cout << "\tNode Index: " << std::distance(runningVectorIndices.begin(), nodeIndex) << std::endl;
		}
		else {
		    curE = *((*distMatrix)[std::distance(runningVectorIndices.begin(), nodeIndex)].rbegin());
		}


	}else {
		curE = (*distMatrix)[node->index][indexCounter];
	}


	curE = curE > maxE ? curE : maxE;

	//std::cout << "Got curE" << std::endl;

	//Check if the node needs inserted at this level
	if(curE < maxEpsilon){
		simp.insert(node->index);
		//Get the largest weight of this simplex
		maxE = curE > node->weight ? curE : node->weight;
		simplexNode* insNode = new simplexNode(simp, maxE);
		insNode->index = curIndex;
		nodeCount++;

		//if depth (i.e. 1 for first iteration) is LT weightGraphSize (starts at 1)
		if(simplexList.size() < simp.size()){
			simplexList.push_back({insNode});
		} else {
			simplexList[simp.size() - 1].insert(insNode);
		}

		//Check if the node has children already...
		if(node->child == nullptr){
			node->child = insNode;
			insNode->parent = node;
			node->children.insert(insNode);

		//Node has children, add to the end of children
		} else {
			//Move to last sibling
			insNode->parent = node;
			insNode->sibling = node->child;
			node->child = insNode;
			node->children.insert(insNode);

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


void simplexTree::printTree(simplexNode* headPointer){
	std::cout << "_____________________________________" << std::endl;
	if(root->child == nullptr){
		std::cout << "Empty tree... " << std::endl;
		return;
	}

	std::cout << "ROOT: " << headPointer->index << "\t" << headPointer << "\t" << headPointer->child << "\t" << headPointer->sibling << std::endl;
	std::cout << "sList Size: " << simplexList.size() << std::endl;

	simplexNode* current;

	for(int i = 0; i < simplexList.size() ; i++){
		std::cout << std::endl << "Dim: " << i << "\tSize: " << simplexList[i].size() << std::endl;
		std::cout << "[index , address, sibling, child, parent]" << std::endl << std::endl;

		for(auto simplexIter = simplexList[i].begin(); simplexIter != simplexList[i].end(); simplexIter++){
			std::cout << (*simplexIter)->index << "\t" << (*simplexIter) << "\t" << (*simplexIter)->sibling << "\t" << (*simplexIter)->child << "\t" << (*simplexIter)->parent << "\t";
			ut.print1DVector((*simplexIter)->simplex);
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

		//Delete the oldest point in the window
		deleteIterative(runningVectorIndices[0]);
		runningVectorIndices.erase(runningVectorIndices.begin() + 0);

		//Create distance matrix row of current vector to each point in the window
		// std::vector<double> distsCurrVec = ut.nearestNeighbors(currentVector, window);
		std::vector<double> distMatrixRow = ut.nearestNeighbors(currentVector, window);

//		//Insert the new point into the distance matrix and complex
//		for(int i = 0; i < (*distMatrix).size(); i++) {
//            (*distMatrix)[i].push_back(distsCurrVec[i+1]);
//		}
//
//		std::vector<double> distMatLastRow( window.size() );
//		distMatrix->push_back(distMatLastRow);

        distMatrix->push_back(distMatrixRow);

		insert();

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

		std::cout << "indexToBeDeleted = " << indexToBeDeleted << '\n';
		deleteIndexRecurse(keyToBeDeleted);
		runningVectorIndices.erase(runningVectorIndices.begin() + indexToBeDeleted);

		insert();

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

//        if(index < 0)  // Impossible event - no deletion
//        {
//            //Delete Row[index]
//            distMatrix->erase(distMatrix->begin() + index);
//
//            //Delete column[index] (row[][index])
//            for(int i = 0; i < (*distMatrix).size(); i++)
//            {
//                if((*distMatrix)[i].size() >= index)
//                    (*distMatrix)[i].erase((*distMatrix)[i].begin() + index);
//            }
//        }

        auto curNodeCount = nodeCount;

		//Delete all entries in the simplex tree
		//printTree(root);
		// std::cout << "simplex->index " << simplex->index << '\n';
		deleteIndexRecurse(simplex->index);

		//printTree(root);

	} else {
		ut.writeDebug("simplexTree","Failed to find vector by index");
	}
	return;
}


void simplexTree::deleteIndexRecurse(int vectorIndex) {
    std::cout << "deleteIndexRecurse vectorIndex = " << vectorIndex << '\n';
    deleteIndexRecurse(vectorIndex, root->child);
    return;
}


void simplexTree::deleteIndexRecurse(int vectorIndex, simplexNode* curNode){
	if(curNode == nullptr){
		std::cout << "Empty tree" << std::endl;
		return;
	}

	//Handle siblings - either they need to be removed (and 'hopped') or recursed
	//	Only need to recurse if the vector could be a member of tree (i.e. < vectorIndex)
	if(curNode->sibling != nullptr && curNode->sibling->index == vectorIndex){

		//Map the current node's sibling to the node following the node for deletion
		simplexNode* tempNode = curNode->sibling;
		curNode->sibling = curNode->sibling->sibling;

		//Delete the orphaned node now that sibling remapping is handled
		//	The sibling still points to the following node, so this is valid
		deleteIndexRecurse(vectorIndex, tempNode);

	} else if(curNode->sibling != nullptr){
		// std::cout << "I am there " << '\n';
		//Recurse to the next sibling and check for deletion
		deleteIndexRecurse(vectorIndex, curNode->sibling);
	}

	//Handle self; if this is the vector index delete branches/nodes
	//	Assume all children / sibling pointers have been alleviated in calling function
	if(curNode->index == vectorIndex){

		//Ensure we aren't the first node of the tree
		if(curNode == root->child){
			root->child = curNode->sibling;
		}

		//Remove our index from the parent
		curNode->parent->children.erase(curNode);

		//Delete this node and sub-nodes in the tree
		deletion(curNode);

	//If this isn't the vectorIndex, need to look at children and remove or recurse
	//	Only need to recurse if the vector could be a member of tree (i.e. < vectorIndex)
	} else if (curNode->child != nullptr && curNode->child->index == vectorIndex){
		simplexNode* tempNode = curNode->child;
		curNode->child = curNode->child->sibling;

		for(auto d : simplexList){
			if((*d.begin()) == tempNode) {
			    std::cout << "d type: " << typeid(d).name() << '\n';

                d.insert(d.begin(), tempNode->sibling);
			}
		}


//        for(auto i = 0; i < simplexList.size(); i++)
//        {
//            if(*(simplexList[i].begin()) == tempNode)
//            {
//                std::cout << "Iteration inside if " << i << '\n';
//
//                simplexList[i].insert(simplexList[i].begin(), tempNode->sibling);
//            }
//        }

		deleteIndexRecurse(vectorIndex, tempNode);

	} else if(curNode->child != nullptr && curNode->child->index < vectorIndex){
	    // std::cout << "I am here again" << '\n';
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

	//Create our new node to insert
	simplexNode* insNode = new simplexNode({indexCounter}, 0);
	insNode->index = indexCounter;

	//Track this index in our current window (for sliding window)
	runningVectorIndices.push_back(insNode);

	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(root == nullptr){
		root = new simplexNode();
		insNode->parent = root;
		root->child = insNode;
		root->children.insert(insNode);
		indexCounter++;
		runningVectorCount++;
		nodeCount++;
		simplexList.push_back({insNode});

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
		recurseInsert((*simplexListIter), indexCounter, 0, 0, {indexCounter});
	}

	//Insert into the right of the tree
	//Child points to last child, and sibling points backwards
	insNode->parent = root;
	insNode->sibling = root->child;
	root->child = insNode;
	root->children.insert(insNode);

	simplexList[0].insert(insNode);

	runningVectorCount++;
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
simplexNode* simplexTree::find(std::set<unsigned>::iterator it, std::set<unsigned>::iterator end, simplexNode* curNode){
	simplexNode* temp = new simplexNode();
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

std::vector<simplexNode*> simplexTree::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<simplexNode*, simplexNode*>& pivotPairs, bool checkEmergent){
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

	//Iterate each child and delete
	auto it = curNode->children.begin();
	while(it != curNode->children.end()){
		deletion(*(it));
		it++;
	}
	nodeCount--;

	//Remove from the simplex list
	simplexList[curNode->simplex.size() - 1].erase(simplexList[curNode->simplex.size() - 1].find(curNode));

	delete curNode;
	return false;
}

void simplexTree::clear(){
	//Clear the simplexTree structure

	if(root != nullptr){
		//Iterate each simplexList[0] entry and delete
		for(auto simp : root->children)
			deletion(simp);

		root = nullptr;
	}

	simplexList.clear();

	simplexOffset = runningVectorCount;
	runningVectorIndices.clear();
	runningVectorCount = 0;
	indexCounter = 0;
	nodeCount = 0;

	return;

}

