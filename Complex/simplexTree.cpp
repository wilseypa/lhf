#include <string>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <iostream>
#include "simplexTree.hpp"

simplexTree::simplexTree(double _maxEpsilon, std::vector<std::vector<double>> _distMatrix, int _maxDim){
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
void simplexTree::recurseInsert(treeNode* node, unsigned curIndex, int depth, double maxE, std::set<unsigned> simp){
	//Incremental insertion
	//Recurse to each child (which we'll use the parent pointer for...)
	treeNode* temp;
	
	isSorted = false;
	double curE = 0;

	//std::cout << runningVectorIndices.size() << "\t" << runningVectorCount << "\t" << indexCounter << "\t" << distMatrix.size() << "\t" << node->index << std::endl;
	if(runningVectorIndices.size() < runningVectorCount+1){
		int offset = runningVectorCount+1 - runningVectorIndices.size();
		if(((int)node->index - offset) > distMatrix.size() || (indexCounter-offset) > distMatrix[(int)node->index - offset].size()){
			std::cout << "DistMatrix access error:" << std::endl;
			std::cout << "DistMatrix size: " << distMatrix.size() << "\tAccess Index: " << ((int)node->index - offset) << std::endl;
			std::cout << "Node Index: " << node->index << "\tOffset: " << offset << std::endl;
			std::cout << "nodeCount: " << nodeCount << "\tindexCount: " << indexCounter << std::endl;
		}
		else
			curE = distMatrix[(int)node->index - offset][indexCounter - offset];

	}else{
		curE = distMatrix[node->index][indexCounter];
	}

	curE = curE > maxE ? curE : maxE;

	//Check if the node needs inserted at this level
	if(curE < maxEpsilon){
		treeNode* insNode = new treeNode();
		insNode->index = curIndex;
		insNode->simplex = simp;
		insNode->simplex.insert(curIndex);
		insNode->sibling = nullptr;
		insNode->child = nullptr;
		nodeCount++;
		simp.insert(node->index);

		//Get the largest weight of this simplex
		insNode->weight = curE > node->weight ? curE : node->weight;
		
		//if depth (i.e. 1 for first iteration) is LT weightGraphSize (starts at 1)
		if(weightEdgeGraph.size() < simp.size()){
			std::vector<std::pair<std::set<unsigned>,double>> tempWEG;
			tempWEG.push_back(std::make_pair(simp, maxE));
			weightEdgeGraph.push_back(tempWEG);
		} else {
			weightEdgeGraph[simp.size() - 1].push_back(std::make_pair(simp, maxE));
		}

		//Check if the node has children already...
		if(node->child == nullptr){
			node->child = insNode;
			insNode->parent = node;

		} else {
			temp = node->child;
			while(temp->sibling != nullptr)
				temp = temp->sibling;
			
			temp->sibling = insNode;
			insNode->parent = temp->parent;
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


void simplexTree::printTree(treeNode* head){
	std::cout << "_____________________________________" << std::endl;
	if(head == nullptr){
		std::cout << "Empty tree... " << std::endl;
		return;
	}

	std::cout << "HEAD: " << head->index << "\t" << head << "\t" << head->child << "\t" << head->sibling << std::endl;

	treeNode* current;
	for(int i = 0; i < dimensions.size() ; i++){
		std::cout << std::endl << "Dim: " << i << std::endl << std::endl;
		current = dimensions[i];

		do{
			std::cout << current->index << "\t" << current << "\t" << current->sibling << "\t" << current->child << "\t" << current->parent << std::endl;
		} while(current->sibling != nullptr && (current = current->sibling) != nullptr);

		std::cout << std::endl;
	}

	std::cout << "_____________________________________" << std::endl;
	return;
}

void simplexTree::insertInductive(){
	//Create our new node to insert
	treeNode* curNode = new treeNode;
	curNode->index = indexCounter;

	//Loop each dimensional list; start at dn, then dn-1, ..., d0 (OLD WAY, INDUCTIVE VR)
	for(int i = (dimensions.size()-1 < maxDimension ? dimensions.size() - 1 : maxDimension); i >= 0; i--){

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
					insNode->sibling = nullptr;

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

	return;
}



//**													**//
//** 				Public Functions 					**//
//**													**//


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

		distMatrix.push_back(distMatrixRow);
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
void simplexTree::deleteIterative(int vectorIndex){
	//Find what row/column of our distance matrix pertain to the vector index
	std::vector<int>::iterator it;
	if((it = std::find(runningVectorIndices.begin(), runningVectorIndices.end(), vectorIndex)) != runningVectorIndices.end()){

		int index = it - runningVectorIndices.begin();

		//Delete the row and column from the distance matrix based on vector index
		if(index > 0){
			
			distMatrix.erase(distMatrix.begin() + index);
		
			for(auto z : distMatrix){
				if(z.size() >= index)
					z.erase(z.begin() + index);
			}
		}

		auto curNodeCount = nodeCount;

		//Delete all entries in the simplex tree with the index...
		deleteIndexRecurse(vectorIndex, head);
		
		deleteWeightEdgeGraph(vectorIndex);
		
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


void simplexTree::deleteIndexRecurse(int vectorIndex, treeNode* curNode){
	//std::cout << "deleteIndexRecurse: " << vectorIndex << std::endl;
	if(curNode == nullptr){
		std::cout << "Empty tree" << std::endl;
		return;
	}
	
	//Handle siblings - either they need to be removed (and 'hopped') or recursed
	if(curNode->sibling != nullptr && curNode->sibling->index == vectorIndex){
		treeNode* tempNode = curNode->sibling;
		curNode->sibling = curNode->sibling->sibling;
		deleteIndexRecurse(vectorIndex, curNode->sibling);

	} else if(curNode->sibling != nullptr){
		deleteIndexRecurse(vectorIndex, curNode->sibling);
	}

	if(curNode->index == vectorIndex){

		if(curNode == head){
			head = curNode->sibling;
			dimensions[0] = curNode->sibling; 
		}
		
		for(auto d : dimensions){
			if(d == curNode)
				d = curNode->sibling;
		}
			
		deletion(curNode);
	} else if (curNode->child != nullptr && curNode->child->index == vectorIndex){
		treeNode* tempNode = curNode->child;
		curNode->child = curNode->child->sibling;
		
		for(auto d : dimensions){
			if(d == tempNode)
				d = tempNode->sibling;
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
	
	if(distMatrix.size() == 0){
		ut.writeDebug("simplexTree","Distance matrix is empty, skipping insertion");
		return;
	}

	//Create our new node to insert
	treeNode* curNode = new treeNode;
	curNode->index = indexCounter;
	std::set<unsigned> tempSet = {curNode->index};
	runningVectorIndices.push_back(indexCounter);

	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(head == nullptr){
		head = curNode;
		indexCounter++;
		runningVectorCount++;
		nodeCount++;
		dimensions.push_back(head);

		std::vector<std::pair<std::set<unsigned>,double>> tempWEG;
		tempWEG.push_back(std::make_pair(tempSet, 0));
		weightEdgeGraph.push_back(tempWEG);

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

	treeNode* temp = dimensions[0];

	//Now let's do it the new way - incremental VR;
	//	Start at the first dimension (d0);
	//		if e (d0, dn) < eps
	//			recurse down tree
	//				recurse
	//			insert to current
	//	iterate to d0->sibling

	do{
		//std::cout << "testDim: " << temp->index << "\t" << dimensions.size() << std::endl;
		recurseInsert(temp, indexCounter, 0, 0, tempSet);
	}while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);

	//Insert into the right of the tree
	temp = dimensions[0];
	treeNode* ins = new treeNode();
	while(temp->sibling != nullptr && (temp = temp->sibling) != nullptr);
	ins->index = indexCounter;
	ins->sibling = nullptr;
	ins->parent = nullptr;
	temp->sibling = ins;
	weightEdgeGraph[0].push_back(std::make_pair(tempSet, 0));

	nodeCount++;
	indexCounter++;
	runningVectorCount++;
	return;
}

void simplexTree::deleteWeightEdgeGraph(int index){
	
	for(unsigned dim = 0; dim < weightEdgeGraph.size(); dim++){
		std::vector<unsigned> d_indices;
		
		for(int vInd = weightEdgeGraph[dim].size()-1; vInd >= 0; vInd--){
			//std::cout << "vInd: " << vInd << std::endl;
			
			if(weightEdgeGraph[dim][vInd].first.find(index) != weightEdgeGraph[dim][vInd].first.end())
				d_indices.push_back(vInd);
		}
		
		int rem = 0;
		for(auto del : d_indices){
			rem++;
			weightEdgeGraph[dim].erase(weightEdgeGraph[dim].begin() + del);
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
	return nodeCount * sizeof(treeNode);
}


std::vector<std::vector<std::pair<std::set<unsigned>,double>>> simplexTree::getAllEdges(double epsilon){

	if(!isSorted){
		for (int i = 0; i < weightEdgeGraph.size(); i++)
			std::sort(weightEdgeGraph[i].begin(), weightEdgeGraph[i].end(), ut.sortBySecond);
		isSorted = true;
	}

	return weightEdgeGraph;
}

void simplexTree::reduceComplex(){
	if(weightEdgeGraph.size() == 0){
		ut.writeDebug("simplexTree","Complex is empty, skipping reduction");
		return;
	}

	//Start with the largest dimension
	ut.writeDebug("simplexTree","Reducing complex, starting simplex count: " + std::to_string(simplexCount()));
	treeNode* cur;

	if(dimensions.size() > 0){
		for(auto i = dimensions.size()-1; i > 1; i--){

			std::vector<std::set<unsigned>> removals;
			std::vector<std::set<unsigned>> checked;

			while(checked.size() != weightEdgeGraph[i].size()){
				cur = dimensions[i];
				do {
					if(std::find(checked.begin(),checked.end(),cur->simplex) == checked.end()){
						auto ret = recurseReduce(std::make_pair(cur->simplex, cur->weight), removals, checked);
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

std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> simplexTree::recurseReduce(std::pair<std::set<unsigned>,double> simplex, std::vector<std::set<unsigned>> removals, std::vector<std::set<unsigned>> checked){
	checked.push_back(simplex.first);
	auto subsets = ut.getSubsets(simplex.first);
	std::set<unsigned> maxFace;

	bool canRemove = true;

	//Look at each face
	for(auto face : subsets){

		//Check if the face is shared; if so, recurse
		for(auto simp : weightEdgeGraph[simplex.first.size() - 1]){

			if(simp != simplex && std::find(checked.begin(), checked.end(), simp.first) == checked.end()){
				auto sDiff = ut.symmetricDiff(simp.first, face,true);

				//This point intersects;
				if (sDiff.size() == 1){
					auto ret = recurseReduce(simp, removals, checked);
					removals = ret.first;
					checked = ret.second;

					//Check if the simplex was not removed
					if(std::find(removals.begin(), removals.end(), simp.first) == removals.end()){
						canRemove = false;
						break;
					}
				}

			}
		}

		//Check if the face is the max face
		double wt = -1;
		if((wt = findWeight(face)) == simplex.second){
			maxFace = face;
		}

	}

	if(canRemove){
		removals.push_back(simplex.first);
		removals.push_back(maxFace);
	}

	return std::make_pair(removals, checked);

}



std::vector<std::vector<unsigned>> simplexTree::getDimEdges(int d,double){
	std::vector<std::vector<unsigned>> ret;

	if(weightEdgeGraph.size() >= d){
		for(auto z : weightEdgeGraph[d]){
			std::vector<unsigned> temp(z.first.begin(), z.first.end());
			ret.push_back(temp);
		}
	}

	return ret;

}

bool simplexTree::find(std::set<unsigned>){

	ut.writeLog("simplexTree","find(std::set<unsigned>) not implemented!");
	return 0;
}

// A recursive function to delete a simplex (and sub-branches) from the tree.
bool simplexTree::deletion(std::set<unsigned> removalEntry) {
	//Remove the entry in the simplex tree
	bool found = true;
	treeNode* curNode = dimensions[0];
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
bool simplexTree::deletion(treeNode* removalEntry) {
	treeNode* curNode = removalEntry;

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

double simplexTree::findWeight(std::set<unsigned> simplex){

	ut.writeLog("simplexTree","findWeight(std::set<unsigned>) not implemented!");
	return 1;
}


std::vector<std::pair<double, std::vector<unsigned>>> simplexTree::getd0Pairs(){
	std::vector<std::pair<double, std::vector<unsigned>>> ret;
	ut.writeLog("simplexTree","getd0Pairs() not implemented!");
	return ret;
}


void simplexTree::clear(){

	//Clear the simplexTree structure
	deletion(head);
	head = nullptr;

	//Clear the weighed edge graph
	for(auto z : weightEdgeGraph){
		z.clear();
	}
	weightEdgeGraph.clear();
	dimensions.clear();


	//runningVectorCount = 0;
	runningVectorIndices.clear();
	//indexCounter = 0;

}

