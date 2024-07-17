#include <string>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <iostream>
// #include <typeinfo>
#include "simplexTree.hpp"
#include <math.h>

template<typename nodeType>
simplexTree<nodeType>::simplexTree(double _maxEpsilon, int _maxDim){
	 
    this->indexCounter = 0;
	this->maxDimension = _maxDim;
	this->maxEpsilon = _maxEpsilon;
	this->simplexType = "simplexTree";
	return;
}

template<typename nodeType>
void simplexTree<nodeType>::outputComplex(){
	return printTree(root);
}

template<typename nodeType>
void simplexTree<nodeType>::recurseInsertDsimplex(simplexTreeNode_P node, std::vector<int> simp,std::vector<std::vector<double>> inputData){
	
	//This algorithm insert a simplex and all its subfaces in the simplex tree. Let σ be a simplex we want to insert with all its subfaces.
	// Let [l0, · · · , lj ] be its word representation. For i from 0 to j we insert, if not already present, a node Nli , storing label li, as a child of the root.
	// We recursively call the algorithm on the subtree rooted at Nli for the insertion of the suffix [li+1,··· ,lj]. 
	// Since the number of subfaces ofr
	//  a simplex of dimension j is 􏰆 􏰁j+1􏰂 = 2j+1, this algorithm takes time O(2jDm).
	sort(simp.begin(),simp.end());
	int firstind =0;
	int lastind = simp.size();
	for(auto x : simp){
        	firstind++;
		std::vector<int> :: const_iterator first = simp.begin() + firstind;
		std::vector<int> :: const_iterator last = simp.begin() + lastind;
		std::vector<int> subsimplex(first,last);
		std::set<unsigned> simplex;
		if(node==nullptr)
	             simplex = {};
		else
		     simplex = node->simpNode->simplex;
		simplex.insert(x);
	    double weight = 0;
		double circumRadius;
		double volume;
		std::vector<double> circumCenter;
		if(simplex.size()>2){
			circumRadius = utils::circumRadius(simplex,this->distMatrix);
			volume = utils::simplexVolume(simplex,this->distMatrix,inputData[0].size());
		}
		else{
			circumRadius = weight/2;
			volume = weight;
		}
		std::cout<<simplex.size()<<std::endl;
		if(simplex.size()>2)
			circumCenter = utils::circumCenter(simplex,inputData);
			
		else if(simplex.size()==2){
 			auto first = simplex.begin();
			std::vector<double> R;
			std::vector<double> A = inputData[*first];
                        std::advance(first, 1);
			std::vector<double> B = inputData[*first];
   	        	std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
		        circumCenter = R;
       }else
   	   circumCenter = inputData[*(simplex.begin())];
 		simplexTreeNode_P insNode = std::make_shared<simplexTreeNode<nodeType>>(simplex, circumRadius);
    		//TODO: dependent on alpha node
            //insNode->simpNode->circumCenter = circumCenter;	
    		//insNode->simpNode->circumRadius = circumRadius;	
       		insNode->simpNode->index = x;
		insNode->simpNode->hash = this->nodeCount;
		this->nodeCount++;
     		if(root == nullptr){
			root = std::make_shared<simplexTreeNode<nodeType>>();
			insNode->parent = root.get();
			root->child = insNode;
			if(subsimplex.size() > 0)
				recurseInsertDsimplex(root->child,subsimplex,inputData);
      	 	}else{
			bool found = false;
			if(node==nullptr){
				node = root;
				insNode->parent = node.get();
				insNode->sibling = node->child;
				node->child = insNode;
				if(subsimplex.size() > 0)
					recurseInsertDsimplex(node->child,subsimplex,inputData);
				}else{
					for(auto it = node->child; it != nullptr; it = it->sibling){
						if(it->simpNode->simplex == simplex){
		    					 found = true;
        				 		 if(subsimplex.size() > 0)
	             				       		  recurseInsertDsimplex(it,subsimplex,inputData);
					  		 break;
       						}
					}
					if(!found){	
	                			if(node->child == nullptr){
	                				node->child = insNode;
	                				insNode->parent = node.get();
                				}else{
							insNode->parent = node.get();
							insNode->sibling = node->child;
       							node->child = insNode;
						}
	           				if(subsimplex.size() > 0)
		        				recurseInsertDsimplex(node->child,subsimplex,inputData);
					}
				}
		}
	}
	return;
}

template<typename nodeType>
void simplexTree<nodeType>::recurseInsert(simplexTreeNode<nodeType>* node, unsigned curIndex, int depth, double maxE, std::set<unsigned> simp){
	//Incremental insertion
	//Recurse to each child (which we'll use the parent pointer for...)
	simplexTreeNode<nodeType>* temp;
	double curE = 0;


	if(this->runningVectorIndices.size() < this->runningVectorCount){

		//Get the positions of the vector in the runningVectorIndices array
		auto nodeIndex = std::find(this->runningVectorIndices.begin(), this->runningVectorIndices.end(), node->simpNode->index);

		if((std::distance(this->runningVectorIndices.begin(), nodeIndex)) > this->distMatrix->size() || (this->indexCounter - (this->runningVectorCount - 1)) > (*this->distMatrix)[std::distance(this->runningVectorIndices.begin(), nodeIndex)].size()){
			std::cout << "DistMatrix access error:" << std::endl;
			std::cout << "\tAttempting to access distMatrix indexes: " << node->simpNode->index << " x " << this->indexCounter << std::endl;
			std::cout << "\tDistMatrix size: " << (*this->distMatrix).size() << std::endl;
			std::cout << "\trviCount: " << this->runningVectorCount << "\t rviSize: " << this->runningVectorIndices.size() << "\tOffset: " << this->simplexOffset << "\tIC: " << this->indexCounter << std::endl;
			std::cout << "\tOffset Indices: " << node->simpNode->index - (this->runningVectorCount - 1) << " x " << this->indexCounter - (this->runningVectorCount - 1) << std::endl;
			std::cout << "\tBackwards size: " << this->distMatrix[this->indexCounter - (this->runningVectorCount - 1)].size() << std::endl;
			std::cout << "\tRow Size: " << this->distMatrix[this->indexCounter - (this->runningVectorCount - 1)].size()  << "\tCurIndex: " << curIndex << std::endl;
			std::cout << "\tNode Index: " << std::distance(this->runningVectorIndices.begin(), nodeIndex) << std::endl;
		}
		else {
		    curE = *((*this->distMatrix)[std::distance(this->runningVectorIndices.begin(), nodeIndex)].rbegin());
		}


	} else {
		curE = (*this->distMatrix)[node->simpNode->index][this->indexCounter];
	}


	curE = curE > maxE ? curE : maxE;

	//Check if the node needs inserted at this level
	if(curE <= this->maxEpsilon){
		simp.insert(node->simpNode->index);
		//Get the largest weight of this simplex
		maxE = curE > node->simpNode->weight ? curE : node->simpNode->weight;
		simplexTreeNode_P insNode = std::make_shared<simplexTreeNode<nodeType>>(simp, maxE);
		insNode->simpNode->index = curIndex;
	        insNode->simpNode->hash = this->nodeCount;
		this->nodeCount++;


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

			temp = insNode->sibling.get();

			//Have to check the children now...
			if(simp.size() <= this->maxDimension){
				do {
					recurseInsert(temp, curIndex, depth + 1, maxE, simp);
				} while((temp = temp->sibling.get()) != nullptr);
			}
		}
	}

	return;
}

template<typename nodeType>
void simplexTree<nodeType>::printTree1(simplexTreeNode_P headPointer){
    
    //Needs to be updated to only alpha nodes, probably specialize printTree()
    
    /*
        if(headPointer == nullptr)
		return;
	if(headPointer->valid){
	for(auto x :headPointer->simpNode->simplex)
		std::cout<<x<<",";
	std::cout<<headPointer->simpNode->filterationvalue<<std::endl;
	}
    	for(auto it = headPointer->child;it!=nullptr;it=it->sibling)
	     printTree1(it);
         */

	return;
}

template<typename nodeType>
void simplexTree<nodeType>:: validateNodes(simplexTreeNode_P headPointer){
//Commenting out temporarily - alphaFilterationValue not part of arguments?
    
     /*   if(headPointer == nullptr)
		return;
	if(headPointer->simpNode->filterationvalue > alphaFilterationValue){
		headPointer->valid = false;
//		std::cout<<headPointer->simpNode->filterationvalue<<"  false\n";
	}
	else{
		headPointer->valid = true;
	
//		std::cout<<headPointer->simpNode->filterationvalue<<" true \n";
	}
//	std::cout<<headPointer->valid<< " "<<alphaFilterationValue;
        headPointer->simpNode->weight = headPointer->simpNode->filterationvalue; 
	for(auto it = headPointer->child;it!=nullptr;it=it->sibling)
	     validateNodes(it);
*/
	return;
}

template<typename nodeType>
void simplexTree<nodeType>::printTree(simplexTreeNode_P headPointer){
	std::cout << "_____________________________________" << std::endl;
	if(root->child == nullptr){
		std::cout << "Empty tree... " << std::endl;
		return;
	}
	auto temp = headPointer.get();

	std::cout << "ROOT: " << headPointer->simpNode->index << "\t" << headPointer << "\t" << headPointer->child << "\t" << headPointer->sibling << std::endl;
	

	std::cout << "[index , address, sibling, child, parent]" << std::endl << std::endl;

	for(auto simplexIter = headPointer; simplexIter != nullptr; simplexIter = simplexIter->sibling){
		std::cout << simplexIter->simpNode->index << "\t";
		std::cout << simplexIter << "\t" ;
		std::cout << simplexIter->sibling << "\t" ;
		std::cout << simplexIter->child << "\t" ;
		std::cout << simplexIter->parent << "\t";
		this->ut.print1DVector(simplexIter->simpNode->simplex);
		
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
		this->ut.print1DVector(simplexIter->simpNode->simplex);
	}
	
	return;
}

// Insert a node into the tree using the distance matrix and a vector index to track changes
template<typename nodeType>
bool simplexTree<nodeType>::insertIterative(std::vector<double> &currentVector, std::vector<std::vector<double>> &window){
	if(window.size() == 0){
		return true;
	}

	if(this->streamEval(currentVector, window)) {   // Point is deemed 'significant'

		//Delete the oldest point in the window
		deleteIterative(this->runningVectorIndices[0]);
		this->runningVectorIndices.erase(this->runningVectorIndices.begin());

		//Create distance matrix row of current vector to each point in the window
		std::vector<double> distsCurrVec = this->ut.nearestNeighbors(currentVector, window);
		// std::vector<double> distMatrixRow = ut.nearestNeighbors(currentVector, window);

		distsCurrVec.erase(distsCurrVec.begin());

		//Insert the new point into the distance matrix and complex
		for(int i = 0; i < (*this->distMatrix).size(); i++) {
            (*this->distMatrix)[i].push_back(distsCurrVec[i]);
		}

		distsCurrVec.push_back(0);
//
//		std::vector<double> distMatLastRow( window.size() );
		this->distMatrix->push_back(distsCurrVec);

//        distMatrix->push_back(distMatrixRow);

		insert();

		this->removedSimplices++;

		return true;
	}

	return false;
}

// Insert a node into the tree using the distance matrix and a vector index to track changes
template<typename nodeType>
bool simplexTree<nodeType>::insertIterative(std::vector<double> &currentVector, std::vector<std::vector<double>> &window, int &keyToBeDeleted, int &indexToBeDeleted){
	if(window.size() == 0){
		return true;
	}

	if(this->streamEval(currentVector, window)) {   // Point is deemed 'significant'

//	     std::cout << "========================== Simplex tree after insertion ==========================" << '\n';
//		 printTree(root);

		std::cout << "indexToBeDeleted = " << indexToBeDeleted << '\n';
		deleteIndexRecurse(keyToBeDeleted);
		this->runningVectorIndices.erase(this->runningVectorIndices.begin() + indexToBeDeleted);

//		 std::cout << "========================== Simplex tree after deletion ==========================" << '\n';
//		 printTree(root);

		insert();

		this->removedSimplices++;

		return true;
	}

	return false;
}


// Delete a node from the tree and from the distance matrix using a vector index
template<typename nodeType>
void simplexTree<nodeType>::deleteIterative(int simplexIndex){

	//Find what row/column of our distance matrix pertain to the vector index
	std::vector<int>::iterator it;
	if((it = std::find(this->runningVectorIndices.begin(), this->runningVectorIndices.end(), simplexIndex)) != this->runningVectorIndices.end()){

		//Index holds the index in the runningVectorIndices array of the simplexNode pointer
		int index = std::distance(this->runningVectorIndices.begin(), it);
		std::cout << "index = " << index << '\n';

		//Delete the row and column from the distance matrix based on vector index
		//	This corresponds to the index into the runnningVectorIndices array

		//Delete Row[index]
		this->distMatrix->erase(this->distMatrix->begin() + index);

		//Delete column[index] (row[][index])
		for(int i = 0; i < (*this->distMatrix).size(); i++) {
            if((*this->distMatrix)[i].size() >= index)
                (*this->distMatrix)[i].erase((*this->distMatrix)[i].begin() + index);
		}


        auto curNodeCount = this->nodeCount;

		//Delete all entries in the simplex tree
		//printTree(root);
		// std::cout << "simplex->index " << simplex->index << '\n';
		deleteIndexRecurse(simplexIndex);

		//printTree(root);

	} else {
		this->ut.writeDebug("simplexTree","Failed to find vector by index");
	}
	return;
}

template<typename nodeType>
void simplexTree<nodeType>::deleteIndexRecurse(int vectorIndex) {
    std::cout << "deleteIndexRecurse vectorIndex = " << vectorIndex << '\n';

    // Since child index is always higher than parent index, no "top node" between root->child and vectorIndex
    // can contain a subtree that has a node->index =  vectorIndex. Therefore, skip to the top node whose sibling
    // has vectorIndex.
    simplexTreeNode<nodeType>* curNode = root->child.get();

    if(curNode->sibling != nullptr)
    {
        while(curNode->sibling->simpNode->index > vectorIndex)
        {
            curNode = curNode->sibling.get();
        }
    }

    deleteIndexRecurse(vectorIndex, curNode);
    return;
}

template<typename nodeType>
void simplexTree<nodeType>::deleteIndexRecurse(int vectorIndex, simplexTreeNode<nodeType>* curNode){
	if(curNode == nullptr){
		std::cout << "Empty tree" << std::endl;
		return;
	}

	//Handle siblings - either they need to be removed (and 'hopped') or recursed
	//	Only need to recurse if the vector could be a member of tree (i.e. < vectorIndex)
	if(curNode->sibling != nullptr && curNode->sibling->simpNode->index == vectorIndex){   // Current node's sibling is to be deleted

		//Map the current node's sibling to the node following the node for deletion
		simplexTreeNode<nodeType>* tempNode = curNode->sibling.get();
		curNode->sibling = curNode->sibling->sibling;

		//Delete the orphaned node now that sibling remapping is handled
		//	The sibling still points to the following node, so this is valid
		deleteIndexRecurse(vectorIndex, tempNode);

	} else if(curNode->sibling != nullptr){
		//Recurse to the next sibling and check for deletion
		deleteIndexRecurse(vectorIndex, curNode->sibling.get());
	}


	//Handle self; if this is the vector index delete branches/nodes
	//	Assume all children / sibling pointers have been alleviated in calling function
	if(curNode->simpNode->index == vectorIndex){
		//Ensure we aren't the first node of the tree
		if(curNode == root->child.get()){
			root->child = curNode->sibling;
		}

		//Remove our index from the parent
		//curNode->parent->children.erase(curNode);

		//Delete this node and sub-nodes in the tree
		deletion(curNode);

	//If this isn't the vectorIndex, need to look at children and remove or recurse
	//	Only need to recurse if the vector could be a member of tree (i.e. < vectorIndex)
	} else if (curNode->child != nullptr && curNode->child->simpNode->index == vectorIndex){
		simplexTreeNode<nodeType>* tempNode = curNode->child.get();
		curNode->child = curNode->child->sibling;

//		for(auto d : simplexList){
//			if((*d.begin()) == tempNode) {
//                d.insert(d.begin(), tempNode->sibling);
//			}
//		}

		deleteIndexRecurse(vectorIndex, tempNode);

	} else if( curNode->child != nullptr && curNode->child->simpNode->index > vectorIndex ) {
	    // If curNode->child->index < vectorIndex, the subtree rooted at curNode->child cannot contain vectorIndex.
		deleteIndexRecurse(vectorIndex, curNode->child.get());

	}

	return;
}


// Insert a node into the tree
template<typename nodeType>
void simplexTree<nodeType>::insert() {
	if(this->distMatrix->size() == 0){
		this->ut.writeDebug("simplexTree","Distance matrix is empty, skipping insertion");
		return;
	}
    std::cout << "insert" << std::endl;
	
	//Create our new node to insert (Ref Count = 1)
    std::set<unsigned> simp({(unsigned)this->indexCounter});;
	simplexTreeNode_P insNode = std::make_shared<simplexTreeNode<nodeType>>(simp, 0); 	
	insNode->simpNode->index = this->indexCounter;

	//Track this index in our current window (for sliding window)
	this->runningVectorIndices.push_back(insNode->simpNode->index);

	//Check if this is the first node (i.e. head)
	//	If so, initialize the head node
	if(root == nullptr){
		root = std::make_shared<simplexTreeNode<nodeType>>();
		insNode->parent = root.get();
		root->child = insNode;
		this->indexCounter++;
		this->runningVectorCount++;
		this->nodeCount++;

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

	this->runningVectorCount++;
	for(auto it = root->child.get(); it != nullptr; it = it->sibling.get()){
		recurseInsert(it, this->indexCounter, 0, 0, {(unsigned)this->indexCounter});
	}

	//Insert into the right of the tree
	// Child points to last child, and sibling points backwards
	insNode->parent = root.get();
	insNode->sibling = root->child;
	root->child = insNode;

	insNode->simpNode->hash = this->nodeCount;
	this->nodeCount++;
	this->indexCounter++;
}


template<typename nodeType>
void simplexTree<nodeType>::deleteWeightEdgeGraph(int index){
    //Commenting out for now, requires weighted edge graph (not in simplexTree)
	/*for(unsigned dim = 0; dim < simplexList.size(); dim++){

		for(auto simplexListIter = simplexList[dim].begin(); simplexListIter != simplexList[dim].end(); ){

			if((*simplexListIter)->simplex.find(index) != (*simplexListIter)->simplex.end())
				simplexList[dim].erase(*simplexListIter);
			else
				simplexListIter++;
		}
	}*/
	return;
}

template<typename nodeType>
int simplexTree<nodeType>::vertexCount(){
	//Return the number of vertices currently represented in the tree
	if(this->runningVectorIndices.size() < this->runningVectorCount+1)
		return this->runningVectorIndices.size();
	return this->indexCounter;
}

template<typename nodeType>
int simplexTree<nodeType>::simplexCount(){
	//Return the number of simplices currently in the tree
	return this->nodeCount;
}

template<typename nodeType>
double simplexTree<nodeType>::getSize(){
	//Size of node: [int + byte (*) + byte (*)] = 18 Bytes
	return this->nodeCount * sizeof(nodeType);
}

//Search for a simplex from a node in the tree
//  Commenting out for now, template is funked

template<typename nodeType>
simplexTree<nodeType>::simplexTreeNode<nodeType>* simplexTree<nodeType>::find(std::set<unsigned>::iterator begin, std::set<unsigned>::iterator end, std::shared_ptr<simplexTreeNode<nodeType>> curNode){
	auto a = curNode.get();
    return find(begin, end, a);
}

template<typename nodeType>
simplexTree<nodeType>::simplexTreeNode<nodeType>* simplexTree<nodeType>::find(std::set<unsigned>::iterator begin, std::set<unsigned>::iterator end, simplexTreeNode<nodeType>* curNode){
	auto it = begin;
    simplexTreeNode<nodeType>* ret;

	while(it != end){
		simplexTreeNode<nodeType>* ptr;
		for(ptr = curNode->child.get(); ptr != nullptr; ptr = ptr->sibling.get()){
			if(ptr->simpNode->index == (*it)){
				++it;
				ret = ptr;
				break;
			}
		}
		if(ptr == nullptr) return nullptr;
	}
	
	return ret;
}

template<typename nodeType>
std::vector<std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>>> simplexTree<nodeType>::getAllEdges(){
	std::vector<std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>>> ret(this->maxDimension + 1, std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>>());
	if(root != nullptr)
		recurseGetEdges(ret, root, 0, this->maxDimension);
	return ret;
}

template<typename nodeType>
void simplexTree<nodeType>::recurseGetEdges(std::vector<std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>>> &edgeList, std::shared_ptr<simplexTreeNode<nodeType>> current, int depth, int maxDepth){
	for(auto ptr = current->child; ptr != nullptr; ptr = ptr->sibling){
		if(ptr->valid)
			edgeList[depth].insert(ptr->simpNode);
		
		if(ptr->child != nullptr && depth+1 <= maxDepth)
			recurseGetEdges(edgeList, ptr, depth+1, maxDepth);
	
	}
	return;
}

template<typename nodeType>
std::vector<nodeType*> simplexTree<nodeType>::getAllCofacets(std::shared_ptr<nodeType> simp){
	return getAllCofacets(simp, std::unordered_map<long long, templateNode_P>(), false);

}


/*
template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexTree<nodeType>::getAllCofacets(const std::set<unsigned>& simp, double simplexWeight, const std::unordered_map<std::shared_ptr<nodeType>, std::shared_ptr<nodeType>>& pivotPairs, bool checkEmergent){
	// Needs rewrite for using std::set<unsigned>& simp instead of pointers
    
    std::vector<templateNode_P> ret;
	/*simplexTreeNode_P parentNode = find(simp->simplex.begin(), simp->simplex.end(), root);
	if(parentNode == nullptr) return ret; //Simplex isn't in the simplex tree

	simplexTreeNode_P tempNode;
	auto it = simp->simplex.end();

	while(true){
		//Insert all of the children in reverse lexicographic order
		for(auto ptr = parentNode->child; ptr != nullptr; ptr = ptr->sibling){
			if(it == simp->simplex.end()) {
					ret.push_back(ptr->simpNode.get()); //All children of simplex are cofacets
			} else {
				tempNode = find(it, simp->simplex.end(), ptr); //See if cofacet is in the tree
				if(tempNode != nullptr&&ptr->valid){
						ret.push_back(tempNode->simpNode.get());
					
					//If we haven't found an emergent candidate and the weight of the maximal cofacet is equal to the simplex's weight
					//		we have identified an emergent pair; at this point we can break because the interval is born and dies at the
					//		same epsilon
					if(checkEmergent && tempNode->simpNode->weight == simp->weight&& this->simplicialComplex != "alpha"){
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
}*/

template<typename nodeType>
std::vector<nodeType*> simplexTree<nodeType>::getAllFacets(std::shared_ptr<nodeType> simp){
	std::vector<nodeType*> ret;
	simplexTreeNode<nodeType>* parentNode = find(simp->simplex.begin(), simp->simplex.end(), root);
	if(parentNode == nullptr) return ret; //Simplex isn't in the simplex tree

	simplexTreeNode<nodeType>* tempNode;
	auto it = simp->simplex.end();

	while(true){
		--it;
		if(parentNode != root.get()) parentNode = parentNode->parent;
		else break;

		//Insert all of the children in reverse lexicographic order
		tempNode = find(std::next(it), simp->simplex.end(), parentNode);
		if(tempNode != nullptr&&tempNode->valid){
			ret.push_back(tempNode->simpNode.get());
		}
	}
	
	return ret;
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexTree<nodeType>::getAllFacets_P(std::shared_ptr<nodeType> simp){
	std::vector<templateNode_P> ret;
	simplexTreeNode<nodeType>* parentNode = find(simp->simplex.begin(), simp->simplex.end(), root);
	if(parentNode == nullptr) return ret; //Simplex isn't in the simplex tree

	simplexTreeNode<nodeType>* tempNode;
	auto it = simp->simplex.end();

	while(true){
		--it;
		if(parentNode != root.get()) parentNode = parentNode->parent;
		else break;

		//Insert all of the children in reverse lexicographic order
		tempNode = find(std::next(it), simp->simplex.end(), parentNode);
		if(tempNode != nullptr&&tempNode->valid){
			ret.push_back(tempNode->simpNode);
		}
	}
	
	return ret;
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexTree<nodeType>::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<std::shared_ptr<nodeType>, std::shared_ptr<nodeType>>& pivotPairs, bool checkEmergent){
	std::vector<templateNode_P> ret;
	auto parentNode = find(simplex.begin(), simplex.end(), root);
	//auto parentNode = parNode.get();
	if(parentNode == nullptr) return ret; //Simplex isn't in the simplex tree
    
    simplexTreeNode<nodeType>* tempNode;
    std::cout << "simplexTree getAllCofacets" << std::endl;
	auto it = simplex.end();
	while(true){
		//Insert all of the children in reverse lexicographic order
		for(auto ptr = parentNode->child; ptr != nullptr; ptr = ptr->sibling){
			if(it == simplex.end()) {
				ret.push_back(ptr->simpNode); //All children of simplex are cofacets
			} else {
				
                //TODO: Reimplement
                tempNode = find(it, simplex.end(), ptr); //See if cofacet is in the tree
				if(tempNode != nullptr){
						ret.push_back(tempNode->simpNode);
					
					//If we haven't found an emergent candidate and the weight of the maximal cofacet is equal to the simplex's weight
					//		we have identified an emergent pair; at this point we can break because the interval is born and dies at the
					//		same epsilon
					if(checkEmergent && tempNode->simpNode->weight == simplexWeight && this->simplicialComplex!="alpha"){
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

template<typename nodeType>
void simplexTree<nodeType>::reduceComplex(){
    
    /* Commenting out; WEG no longer used
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
    */
	return;
}


template<typename nodeType>
void simplexTree<nodeType>:: buildAlphaComplex(std::vector<std::vector<int>> dsimplexmesh, int npts,std::vector<std::vector<double>> inputData){
    
    //Commenting out for now, need to define cmpByWeight
    /*
	std::set<simplexNode_P,cmpByWeight> dsimplexes;
for(auto simplex : dsimplexmesh){
	std::set<unsigned> simplexset(simplex.begin(),simplex.end());
       	recurseInsertDsimplex(root, simplex,inputData);
	simplexNode_P simp = std::make_shared<simplexNode>(simplexNode(simplexset,0));
		dsimplexes.insert(simp);

}


//validateNodes(root);
//printTree1(root);
//ALPHA COMPLEX FILTERTION BASED on FOLLOWING algorithm
/*Filtration value computation algorithm


for i : dimension →0 do
   for all σ of dimension i
        if filtration(σ) is NaN then
            filtration(σ)=α2(σ)
        end if
        for all τ face of σ do            // propagate alpha filtration value
          if  filtration(τ) is not NaN then
               filtration(τ) = min( filtration(τ), filtration(σ) )
          else
             if τ is not Gabriel for σ then
               filtration(τ) = filtration(σ)
          end if
       end if
     end for
  end for
end for

make_filtration_non_decreasing()
prune_above_filtration()

  

std::vector<std::set<simplexNode_P, cmpByWeight>>  edges = getAllEdges();
for(auto dim = edges.size()-1;dim >0;dim--){
	for(auto simp : edges[dim]){
		simplexTreeNode *simpTN = find(simp->simplex.begin(),simp->simplex.end(),root);
		if(simpTN->simpNode->filterationvalue == -1)
			simpTN->simpNode->filterationvalue = simpTN->simpNode->circumRadius;

	      if(dim>0){
		std::vector<simplexNode_P> facets = getAllFacets_P(simp);
                for(auto face : facets){
			bool gabriel = true;
			std::vector<unsigned> points_check(simp->simplex.size());
			std::vector<unsigned> guilty_points_check;
                        if(face->filterationvalue !=-1){
				face->filterationvalue = std::min(face->filterationvalue ,simpTN->simpNode->filterationvalue);
			}
		        else {
				std::vector<unsigned>::iterator it;
				it=std::set_difference (simp->simplex.begin(), simp->simplex.end(), face->simplex.begin(), face->simplex.end(), points_check.begin());
				points_check.resize(it-points_check.begin());
				for(it = points_check.begin(); it !=points_check.end();++it){
					std::vector<double> coordinates;
					for(int i =0;i<face->circumCenter.size();i++)
						coordinates.push_back(inputData[*it][i]);
					double distance = utils::vectors_distance(coordinates,face->circumCenter);
					if(pow(distance,2)<face->circumRadius){
						gabriel = false;
						guilty_points_check.push_back((*it));
					}
				}
			}
			if(!gabriel){
				std::vector<unsigned> v(simp->simplex.size());
				std::vector<unsigned>::iterator it;
				it=std::set_union (face->simplex.begin(), face->simplex.end(), guilty_points_check.begin(), guilty_points_check.end(), v.begin());
				v.resize(it-v.begin());
				std::set<unsigned> simplexsetNG(v.begin(),v.end());
	                	simplexTreeNode *cofacetNG = find(simplexsetNG.begin(),simplexsetNG.end(),root);
		        	face->filterationvalue = cofacetNG->simpNode->filterationvalue;
			}
		}
	      }
	}
        for(auto simp: edges[0]){
		simplexTreeNode *simpTN = find(simp->simplex.begin(),simp->simplex.end(),root);
		if(simpTN->simpNode->filterationvalue == -1)
			simpTN->simpNode->filterationvalue = simpTN->simpNode->circumRadius;
	}
  
}
//Reinserting to sort by filterationvalue and remove simplexes with weight greater than alphafilteration value
validateNodes(root);*/
return;
}

template<typename nodeType>
std::pair<std::vector<std::set<unsigned>>, std::vector<std::set<unsigned>>> simplexTree<nodeType>::recurseReduce(simplexTreeNode_P simplex, std::vector<std::set<unsigned>> removals, std::vector<std::set<unsigned>> checked){
	checked.push_back(simplex->simpNode->simplex);
	auto subsets = this->ut.getSubsets(simplex->simpNode->simplex);
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

template<typename nodeType>
bool simplexTree<nodeType>::find(std::set<unsigned>){

	this->ut.writeLog("simplexTree","find(std::set<unsigned>) not implemented!");
	return 0;
}

// A recursive function to delete a simplex (and sub-branches) from the tree.
template<typename nodeType>
bool simplexTree<nodeType>::deletion(std::set<unsigned> removalEntry) {
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


	this->ut.writeLog("simplexTree","deletion(std::set<unsigned>) not implemented!");
	return false;

}


// A recursive function to delete a simplex (and sub-branches) from the tree.
template<typename nodeType>
bool simplexTree<nodeType>::deletion(simplexTreeNode<nodeType>* removalEntry) {
	simplexTreeNode<nodeType>* curNode = removalEntry;

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

template<typename nodeType>
void simplexTree<nodeType>::clear(){
	//Clear the simplexTree structure
	root = nullptr;

	this->simplexOffset = this->runningVectorCount;
	this->runningVectorIndices.clear();
	this->runningVectorCount = 0;
	this->indexCounter = 0;
	this->nodeCount = 0;

	return;

}

//Explicit Template Class Instantiation
template class simplexTree<simplexNode>;
template class simplexTree<alphaNode>;
template class simplexTree<witnessNode>;

