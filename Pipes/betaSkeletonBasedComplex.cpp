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
	std::set<std::vector<unsigned>> dsimplexmesh;
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
 //   std::cout<<averageDistance<<" ";
    count =0;
    std::vector<std::pair<int,int>> neighborsepsilon;
    
   
    for(unsigned index = 0; index < inData.inputData.size(); index++){
		    std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); //All neighbors in epsilon-ball
		    int n = neighbors.size();
		   // neighbors.erase(std::remove_if(neighbors.begin(),neighbors.end(),[&index](size_t x){return x<=index;}),neighbors.end());
            neighborsepsilon.push_back(std::make_pair(n,index));
			
    }

   
    sort(neighborsepsilon.begin(),neighborsepsilon.end());
    std::vector<int> toremove;
	for(unsigned indexold = 0; indexold < inData.inputData.size(); indexold++){
		    toremove.push_back(neighborsepsilon[indexold].second);
		    unsigned index = neighborsepsilon[indexold].second;
		    std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); //All neighbors in epsilon-ball
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
				dsimplexIndexed.push_back(difference[dsimplex[i]-1]);

			if(checkInsertDsimplex(dsimplexIndexed,inData,this->beta,averageDistance,tree)){	
					 std::sort(dsimplexIndexed.begin(), dsimplexIndexed.end());		  
				     dsimplexmesh.insert(dsimplexIndexed);
					 count++;
		   }
            
	
	while((*first) != n-dim+1){  	

  	    std::vector<unsigned>::iterator mt = last;
        while (*(--mt) == n-(last-mt)+1);
        (*mt)++;
        while (++mt != last) *mt = *(mt-1)+1;

        std::vector<unsigned> dsimplexIndexed1;

       dsimplexIndexed1.push_back(index);
       for(int i=0;i<dim;i++){
			dsimplexIndexed1.push_back(difference[dsimplex[i]-1]);
			//std::cout<<dsimplexIndexed1[i]<<" \n";
		}
            if(checkInsertDsimplex(dsimplexIndexed1,inData,this->beta,averageDistance,tree)){
			std::sort(dsimplexIndexed1.begin(), dsimplexIndexed1.end());
           	     	dsimplexmesh.insert(dsimplexIndexed1);
		    	count++;	    	
	       
	    		}
    	 	}      
	       
	  }
	}
	
	((alphaComplex<nodeType>*)inData.complex)->buildBetaComplex(dsimplexmesh,inData.inputData.size(),inData.inputData,this->beta,this->betaMode);
	std::ofstream file("PHdSphereDimensionWiseMeshSize.txt",std::ios_base::app);
        file<<this->betaMode<<","<<inData.inputData.size()<<","<<inData.inputData[0].size()<<","<<this->beta<<","<<dsimplexmesh.size()<<std::endl;	
	file.close();
	this->ut.writeDebug("betaSkeletonBasedComplex Pipe", "\tbetaSkeletonBasedComplex Size: ");
	return;
}

template <>

unsigned betaSkeletonBasedComplex<alphaNode>:: selectCenter(std::vector<double> hpcofffaces, std::vector<std::vector<double>> betaCenters,std::vector<double> otherPoint){
   double valuebetaCenter1 =0, valuebetaCenter2=0, valueotherPoint=0;
   for(unsigned i =0;i<hpcofffaces.size();i++){
	   valuebetaCenter1 += hpcofffaces[i]*betaCenters[0][i];
	   valuebetaCenter2 += hpcofffaces[i]*betaCenters[1][i];
	   valueotherPoint += hpcofffaces[i]*otherPoint[i];
	}
	   if(valueotherPoint < 0){
			if(valuebetaCenter1 <0)
				return 0;
			else if(valuebetaCenter2 <0)
				return 1;
		}
		else if(valueotherPoint >0){
			if(valuebetaCenter1 >0)
				return 0;
			else if(valuebetaCenter2 >0)
				return 1;
		}
return 1;
}
template <>

bool betaSkeletonBasedComplex<alphaNode>:: checkInsertDsimplex(std::vector<unsigned> dsimplex,pipePacket<alphaNode> &inData,double beta,double averageDistance,kdTree tree){
	double maxEdge = 0;
	for(auto x : dsimplex)
		for(auto y : dsimplex)
			if(maxEdge < (*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[x][y])
				maxEdge = (*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[x][y];
				
	if(maxEdge > this->epsilon)
	    return false;
	            
	std::vector<size_t> neighborsfinalLune;
        
	bool intersection;
	if(beta <0)
		exit(0);
	else if(beta==0)
		return true;
    	else if(beta <= 1){
		intersection = true;
		beta = 1/beta;
	}
	else
		intersection = false;

			if(this->betaMode == "highDimCircle"){
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
   			double circumRadius;
   			if(simplex.size()>2)
   				circumRadius = utils::circumRadius(simplex,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
   			else
				circumRadius = pow((*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[dsimplex[0]][dsimplex[1]]/2,2);
			std::vector<double> hpcoff = utils::nullSpaceOfMatrix(simplex,inData.inputData,circumCenter,sqrt(circumRadius));
			std::vector<std::vector<double>> betaCenters = utils::betaCentersCalculation(hpcoff, beta, sqrt(circumRadius),circumCenter);
			double betaRadius = beta*sqrt(circumRadius);
		//	double volume = utils::simplexVolume(simplex,((alphaComplex<alphaNode>*)inData.complex)->distMatrix,inData.inputData[0].size());

		    std::vector<size_t> neighbors1 = tree.neighborhoodIndices(betaCenters[0], betaRadius); //All neighbors in epsilon-ball
		    std::vector<unsigned> toremove(simplex.begin(),simplex.end()); 
		    std::sort(toremove.begin(),toremove.end());
            std::sort(neighbors1.begin(), neighbors1.end());
            std::vector<int> neg1;
            std::set_difference(neighbors1.begin(),neighbors1.end(),toremove.begin(),toremove.end(), std::back_inserter(neg1));
		    
    		std::vector<size_t> neighbors2 = tree.neighborhoodIndices(betaCenters[1], betaRadius); //All neighbors in epsilon-ball
            std::sort(neighbors2.begin(), neighbors2.end());
            std::vector<int> neg2;
            std::set_difference(neighbors2.begin(),neighbors2.end(),toremove.begin(),toremove.end(), std::back_inserter(neg2));
    		
    		
    		
            std::vector<size_t> neighbors;
            std::sort (neg1.begin(),neg1.end());
			std::sort (neg2.begin(),neg2.end()); 
			
            if(intersection == true){
				std::vector<size_t> v(std::min(neg1.size(),neg2.size()));
				std::vector<size_t>::iterator it;
				
				it=std::set_intersection (neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
				v.resize(it-v.begin()); 
				neighbors = v;				
			}
			else{
				std::vector<size_t> v(neg1.size() + neg2.size());
				std::vector<size_t>::iterator it;	
		
				it=std::set_union (neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
				v.resize(it-v.begin()); 
				neighbors = v;
			}
		
			if(neighbors.size() == 0){
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
			}
			else if(this->betaMode == "highDimLune"){

			       
				std::set<unsigned> simplex(dsimplex.begin(),dsimplex.end());
			     
			     
				 std::vector<double> circumCenter;
				 std::vector<double> circumCenterfaces;
				 std::vector<double> circumCenterfaces1;

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
				double circumRadius;
				if(simplex.size()>2)
					circumRadius = utils::circumRadius(simplex,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
				else
					circumRadius = pow((*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[dsimplex[0]][dsimplex[1]]/2,2);
				bool first = true;
				for (auto x : simplex){
					std::vector<double> betaCenter;
					for(unsigned y =0 ;y< inData.inputData[0].size();y++)
						betaCenter.push_back(beta*circumCenter[y] + (1-beta)*inData.inputData[x][y]);
					
		           double betaRadius = beta*sqrt(circumRadius);

		    
					std::vector<size_t> neighbors1faces1 = tree.neighborhoodIndices(betaCenter, betaRadius); //All neighbors in epsilon-ball
					neighbors1faces1.erase(std::remove(neighbors1faces1.begin(),neighbors1faces1.end(),x),neighbors1faces1.end());
					if(!first){
						std::sort (neighborsfinalLune.begin(),neighborsfinalLune.end());
						std::sort (neighbors1faces1.begin(),neighbors1faces1.end()); 
			
						std::vector<size_t> v1(std::min(neighborsfinalLune.size(),neighbors1faces1.size()));
						std::vector<size_t>::iterator it1;
				
						it1=std::set_intersection (neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
						v1.resize(it1-v1.begin()); 
						neighborsfinalLune = v1;
					}	
					else
						neighborsfinalLune = neighbors1faces1;
                    first= false;
           
				}
				
	/*			std::vector<double> hpcoff = utils::nullSpaceOfMatrix(simplex,inData.inputData,circumCenter,sqrt(circumRadius));
		
			     int dim = dsimplex.size()-1;
			     int n = dim+1;
				 std::vector<unsigned> dsimplex1(dim);
				 std::vector<unsigned> dsimplexIndexed;

				std::vector<unsigned>::iterator first = dsimplex1.begin(), last = dsimplex1.end();

				std::generate(first, last, UniqueNumber);
			
			for(int i=0;i<dim;i++){
				dsimplexIndexed.push_back(dsimplex[dsimplex1[i]-1]);
			}
                std::set<unsigned> dsimplexIndexedset(dsimplexIndexed.begin(),dsimplexIndexed.end());

				if(dsimplexIndexed.size()>2)
					circumCenterfaces = utils::circumCenter(dsimplexIndexedset,inData.inputData);
				else if(dsimplexIndexed.size()==2){
 				auto first = dsimplexIndexed.begin();
				std::vector<double> R;
				std::vector<double> A = inData.inputData[*first];
      			std::advance(first, 1);
				std::vector<double> B = inData.inputData[*first];
	   			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
				circumCenterfaces = R;
				}
				double circumRadiusfaces;
				if(simplex.size()>2)
					circumRadiusfaces = utils::circumRadius(dsimplexIndexedset,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
				else
					circumRadiusfaces = pow((*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[dsimplexIndexed[0]][dsimplexIndexed[1]]/2,2);
				std::vector<double> hpcofffaces = utils::nullSpaceOfMatrix(dsimplexIndexedset,inData.inputData,circumCenterfaces,sqrt(circumRadiusfaces));
				std::vector<std::vector<double>> betaCentersfaces = utils::betaCentersCalculation(hpcofffaces, beta, sqrt(circumRadiusfaces),circumCenterfaces);
				double betaRadiusfaces;
				if(circumRadiusfaces > circumRadius)
					betaRadiusfaces = beta*sqrt(circumRadiusfaces);
				else
					betaRadiusfaces = beta*sqrt(circumRadius);

				//Decide on beta center
				// extract the other points
		    std::sort(dsimplexIndexed.begin(),dsimplexIndexed.end());
            std::sort(dsimplex.begin(), dsimplex.end());
            std::vector<unsigned> otherPoint;
            std::set_difference(dsimplex.begin(),dsimplex.end(),dsimplexIndexed.begin(),dsimplexIndexed.end(), std::back_inserter(otherPoint));
		    
		    auto selectedCenter = selectCenter(hpcofffaces,betaCentersfaces,inData.inputData[otherPoint[0]]);
		    
		    std::vector<size_t> neighbors1faces = tree.neighborhoodIndices(betaCentersfaces[selectedCenter], betaRadiusfaces); //All neighbors in epsilon-ball
		    std::vector<unsigned> toremovefaces(dsimplexIndexed.begin(),dsimplexIndexed.end()); 
		    std::sort(toremovefaces.begin(),toremovefaces.end());
            std::sort(neighbors1faces.begin(), neighbors1faces.end());
            std::vector<size_t> neg1faces;
            std::set_difference(neighbors1faces.begin(),neighbors1faces.end(),toremovefaces.begin(),toremovefaces.end(), std::back_inserter(neg1faces));
		    neighborsfinalLune = neg1faces;
			
			while((*first) != n-dim+1){  	

			std::vector<unsigned>::iterator mt = last;
			while (*(--mt) == n-(last-mt)+1);
				(*mt)++;
			while (++mt != last) *mt = *(mt-1)+1;

			std::vector<unsigned> dsimplexIndexed1;

			for(int i=0;i<dim;i++){
				dsimplexIndexed1.push_back(dsimplex[dsimplex1[i]-1]);
			}
				std::set<unsigned> dsimplexIndexed1set(dsimplexIndexed1.begin(),dsimplexIndexed1.end());

				if(dsimplexIndexed1.size()>2)
					circumCenterfaces1 = utils::circumCenter(dsimplexIndexed1set,inData.inputData);
				else if(dsimplexIndexed1.size()==2){
 				auto first = dsimplexIndexed1.begin();
				std::vector<double> R;
				std::vector<double> A = inData.inputData[*first];
      			std::advance(first, 1);
				std::vector<double> B = inData.inputData[*first];
	   			std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R),[](double e1,double e2){return ((e1+e2)/2);});
				circumCenterfaces1 = R;
				}
				double circumRadiusfaces1;
				if(simplex.size()>2)
					circumRadiusfaces1 = utils::circumRadius(dsimplexIndexed1set,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
				else
					circumRadiusfaces1 = pow((*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[dsimplexIndexed1[0]][dsimplexIndexed1[1]]/2,2);
				std::vector<double> hpcofffaces1 = utils::nullSpaceOfMatrix(dsimplexIndexed1set,inData.inputData,circumCenterfaces1,sqrt(circumRadiusfaces1));
				std::vector<std::vector<double>> betaCentersfaces1 = utils::betaCentersCalculation(hpcofffaces1, beta, sqrt(circumRadiusfaces1),circumCenterfaces1);
				double betaRadiusfaces1 = beta*sqrt(circumRadiusfaces1);

				//Decide on beta center
				// extract the other points
		    std::sort(dsimplexIndexed1.begin(),dsimplexIndexed1.end());
            std::sort(dsimplex.begin(), dsimplex.end());
            std::vector<unsigned> otherPoint1;
            std::set_difference(dsimplex.begin(),dsimplex.end(),dsimplexIndexed1.begin(),dsimplexIndexed1.end(), std::back_inserter(otherPoint1));
		    
		    auto selectedCenter1 = selectCenter(hpcofffaces1,betaCentersfaces1,inData.inputData[otherPoint1[0]]);
		    
		    std::vector<size_t> neighbors1faces1 = tree.neighborhoodIndices(betaCentersfaces1[selectedCenter1], betaRadiusfaces1); //All neighbors in epsilon-ball
		    std::vector<unsigned> toremovefaces1(dsimplexIndexed1.begin(),dsimplexIndexed1.end()); 
		    std::sort(toremovefaces1.begin(),toremovefaces1.end());
            std::sort(neighbors1faces1.begin(), neighbors1faces1.end());
            std::vector<int> neg1faces1;
            std::set_difference(neighbors1faces1.begin(),neighbors1faces1.end(),toremovefaces1.begin(),toremovefaces1.end(), std::back_inserter(neg1faces1));
		    
            std::sort (neighborsfinalLune.begin(),neighborsfinalLune.end());
			std::sort (neg1faces1.begin(),neg1faces1.end()); 
			
			std::vector<size_t> v1(std::min(neighborsfinalLune.size(),neg1faces1.size()));
			std::vector<size_t>::iterator it1;
				
			it1=std::set_intersection (neg1faces1.begin(), neg1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
			v1.resize(it1-v1.begin()); 
			neighborsfinalLune = v1;				
     	}
     	* */
	//	std::cout<<count<<" "<<dsimplexmesh.size()<<std::endl;
		}
	//	std::cout<<"neighborsfinalLune.size() :: "<<neighborsfinalLune.size()<<"\n";
		if(this->betaMode == "highDimLune" && neighborsfinalLune.size() == 0){
				   //recurseInsertDsimplex(root, dsimplex,inputData);
				return true;
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
	
	pipe = configMap.find("betaMode");
	if(pipe != configMap.end())
		this->betaMode = configMap["betaMode"].c_str();
		
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		this->epsilon = std::atof(configMap["epsilon"].c_str());
		
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
