/*
 * betaSkeletonBasedComplexPipe hpp + cpp extend the basePipe class for calculating the 
 * beta Skeleton Based Complex generation for data input
 * 
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <set>
#include <algorithm>
#include "betaSkeletonBasedComplex.hpp"
#include "alphaComplex.hpp"

#include "utils.hpp"

// basePipe constructor
template <typename nodeType>
betaSkeletonBasedComplex<nodeType>::betaSkeletonBasedComplex(){
	this->pipeType = "betaSkeletonBasedComplex";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void betaSkeletonBasedComplex<nodeType>::runPipe(pipePacket<nodeType> &inData){
	// Generate Beta Skeleton Based Complex
	
	// Temporarily commenting this out - need to check inData.complex type
	//		If type is the graph-based simplexArrayList (inherited) then
	//			cast to gAL and run non-virtual function:
	
	//	((graphArrayList*)inData.complex)->graphInducedComplex(dim,inData.inputData,beta);
	std::vector<std::vector<unsigned>> dsimplexmesh;
	kdTree tree(inData.inputData, inData.inputData.size()); //KDTree for efficient nearest neighbor search
	int dim = inData.inputData[0].size();
    double distanceSum = 0;
	int count = 0;
	//Computing the average point cloud distance
	for(int x =0;x < inData.inputData.size();x++)
		for(int y = x+1; y < inData.inputData.size();y++){
			distanceSum += (*((alphaComplex<nodeType>*)inData.complex)->distMatrix)[x][y];
			count= count+1;
		}
	double averageDistance = distanceSum/count;
    std::cout<<averageDistance<<" ";
    count =0;
    std::vector<std::pair<int,int>> neighborsepsilon;
    
    for(unsigned index = 0; index < inData.inputData.size(); index++){
		    std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->enclosingRadius); //All neighbors in epsilon-ball
		    int n = neighbors.size();
		   // neighbors.erase(std::remove_if(neighbors.begin(),neighbors.end(),[&index](size_t x){return x<=index;}),neighbors.end());
            neighborsepsilon.push_back(std::make_pair(n,index));
			
    }
    sort(neighborsepsilon.begin(),neighborsepsilon.end());
    std::vector<int> toremove;
	for(unsigned indexold = 0; indexold < inData.inputData.size(); indexold++){
		    toremove.push_back(neighborsepsilon[indexold].second);
		    unsigned index = neighborsepsilon[indexold].second;
		    std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->enclosingRadius); //All neighbors in epsilon-ball
		    int n = neighbors.size();
//		    neighbors.erase(std::remove_if(neighbors.begin(),neighbors.end(),[&index](size_t x){return x<=index;}),neighbors.end());
            std::sort(toremove.begin(),toremove.end());
            std::sort(neighbors.begin(), neighbors.end());
            std::vector<int> difference;
            std::set_difference(neighbors.begin(),neighbors.end(),toremove.begin(),toremove.end(), std::back_inserter(difference));
			n = difference.size();
		    
			if(n>=dim){
			std::vector<unsigned> dsimplex(dim);
			std::vector<unsigned> dsimplexIndexed;

			std::vector<unsigned>::iterator first = dsimplex.begin(), last = dsimplex.end();

			std::generate(first, last, UniqueNumber);
			
			dsimplexIndexed.push_back(index);
			for(int i=0;i<dim;i++)
				dsimplexIndexed.push_back(neighbors[dsimplex[i]-1]);

			if(checkInsertDsimplex(dsimplexIndexed,inData,this->beta,averageDistance,tree)){			     
				     dsimplexmesh.push_back(dsimplexIndexed);
					 count++;
		   }
            
	while((*first) != n-dim+1){  	

  	    std::vector<unsigned>::iterator mt = last;
        while (*(--mt) == n-(last-mt)+1);
        (*mt)++;
        while (++mt != last) *mt = *(mt-1)+1;

        std::vector<unsigned> dsimplexIndexed1;

       dsimplexIndexed1.push_back(index);
       for(int i=0;i<dim;i++)
			dsimplexIndexed1.push_back(neighbors[dsimplex[i]-1]);
		
            if(checkInsertDsimplex(dsimplexIndexed1,inData,this->beta,averageDistance,tree)){
                dsimplexmesh.push_back(dsimplexIndexed1);
		    	count++;	
		    	
	    	}
     	}
		std::cout<<count<<std::endl;
	  }
	}
	((alphaComplex<nodeType>*)inData.complex)->buildBetaComplex(dsimplexmesh,inData.inputData.size(),inData.inputData);
	
	this->ut.writeDebug("betaSkeletonBasedComplex Pipe", "\tbetaSkeletonBasedComplex Size: ");
	return;
}
template <>

bool betaSkeletonBasedComplex<alphaNode>:: checkInsertDsimplex(std::vector<unsigned> dsimplex,pipePacket<alphaNode> &inData,double beta,double averageDistance,kdTree tree){
	double maxEdge = 0;
	for(auto x : dsimplex)
		for(auto y : dsimplex)
			if(maxEdge < (*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[x][y])
				maxEdge = (*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[x][y];
				
	if(maxEdge > this->enclosingRadius)
	    return false;
	
	bool intersection;
    if(beta < 1){
		intersection = true;
		beta = 1/beta;
	}
	else
		intersection = false;

   
    std::set<unsigned> simplex(dsimplex.begin(),dsimplex.end());

            std::vector<double> circumCenter;

            if(simplex.size()>2)
				circumCenter = utils::circumCenter(simplex,inData.inputData);
			else if(simplex.size()==2){
 				auto first = simplex.begin();
				std::vector<double> R;
				std::vector<double> A = inData.inputData[*first];
      			std::advance(first, 1);
				std::vector<double> B = inData.inputData[*first];
	   			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
				circumCenter = R;
   			 }
			double circumRadius = utils::circumRadius(simplex,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
			std::vector<double> hpcoff = utils::nullSpaceOfMatrix(simplex,inData.inputData,circumCenter,sqrt(circumRadius));
			std::vector<std::vector<double>> betaCenters = utils::betaCentersCalculation(hpcoff, beta, sqrt(circumRadius),circumCenter);
			double betaRadius = beta*sqrt(circumRadius);
			double volume = utils::simplexVolume(simplex,((alphaComplex<alphaNode>*)inData.complex)->distMatrix,inData.inputData[0].size());

		    std::vector<size_t> neighbors1 = tree.neighborhoodIndices(betaCenters[0], betaRadius); //All neighbors in epsilon-ball
    		std::vector<size_t> neighbors2 = tree.neighborhoodIndices(betaCenters[1], betaRadius); //All neighbors in epsilon-ball
            std::vector<size_t> neighbors;
            if(intersection == true){
				std::vector<size_t> v(std::min(neighbors1.size(),neighbors2.size()));
				std::vector<size_t>::iterator it;
				std::sort (neighbors1.begin(),neighbors1.end());
				std::sort (neighbors2.begin(),neighbors2.end()); 
				it=std::set_intersection (neighbors1.begin(), neighbors1.end(), neighbors2.begin(), neighbors2.end(), v.begin());
				v.resize(it-v.begin()); 
				neighbors = v;				
			}
			else{
				std::vector<size_t> v(neighbors1.size() + neighbors2.size());
				std::vector<size_t>::iterator it;	
				std::sort (neighbors1.begin(),neighbors1.end());
				std::sort (neighbors2.begin(),neighbors2.end()); 
				it=std::set_union (neighbors1.begin(), neighbors1.end(), neighbors2.begin(), neighbors2.end(), v.begin());
				v.resize(it-v.begin()); 
				neighbors = v;
			}
			if(neighbors.size() <= simplex.size()){
				    //recurseInsertDsimplex(root, dsimplex,inputData);
					return true;
/*			    	simplexNode_P tot = std::make_shared<simplexNode>(simplexNode((*it)->simplex, betaRadius));
					tot->simplex.insert(pt);
					tot->hash = (*it)->hash + bin.binom(pt, tot->simplex.size());
					tot->circumCenter = circumCenter;	
					tot->circumRadius = circumRadius;	
					tot->betaCenters = betaCenters;
					tot->betaRadius = betaRadius;
					tot->hpcoff = hpcoff;
					tot->volume = volume;
					simplexList[d].insert(tot);   //insert simplex into complex
*/
				}
return false;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool betaSkeletonBasedComplex<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();
	
	pipe = configMap.find("beta");
	if(pipe != configMap.end())
		this->beta = std::atof(configMap["beta"].c_str());
		
	this->ut = utils(strDebug, this->outputFile);
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		this->enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else return false;

	this->configured = true;
	this->ut.writeDebug("betaSkeletonBasedComplex Pipe ","Configured with parameters { eps: " + configMap["epsilon"] + configMap["beta"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");
	
	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void betaSkeletonBasedComplex<nodeType>::outputData(pipePacket<nodeType> &inData){
	// Output related to betaSkeletonBasedComplex
	return;
}
template <typename nodeType>
bool betaSkeletonBasedComplex<nodeType>:: checkInsertDsimplex(std::vector<unsigned> dsimplex,pipePacket<nodeType> &inData,double beta,double averageDistance,kdTree tree){
	std::cout<<"No function Defined";
	return false;
}

template class betaSkeletonBasedComplex<simplexNode>;
template class betaSkeletonBasedComplex<alphaNode>;
template class betaSkeletonBasedComplex<witnessNode>;
