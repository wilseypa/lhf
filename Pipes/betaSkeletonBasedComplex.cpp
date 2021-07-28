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
     	count =0;
    	std::vector<std::pair<int,int>> neighborsepsilon;

   	for(unsigned index = 0; index < inData.inputData.size(); index++){
    		std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); //All neighbors in epsilon-ball
		int n = neighbors.size();
                neighborsepsilon.push_back(std::make_pair(n,index));
       	}
       	sort(neighborsepsilon.begin(),neighborsepsilon.end());
    	std::vector<int> toremove;
	for(unsigned indexold = 0; indexold < inData.inputData.size(); indexold++){
		    toremove.push_back(neighborsepsilon[indexold].second);
		    unsigned index = neighborsepsilon[indexold].second;
		    std::vector<size_t> neighbors = tree.neighborhoodIndices(inData.inputData[index], this->epsilon); //All neighbors in epsilon-ball
		    int n = neighbors.size();
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
				    	    dsimplexIndexed1.push_back(difference[dsimplex[i]-1]);
				    if(checkInsertDsimplex(dsimplexIndexed1,inData,this->beta,averageDistance,tree)){
					    std::sort(dsimplexIndexed1.begin(), dsimplexIndexed1.end());
           	     			    dsimplexmesh.push_back(dsimplexIndexed1);
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
        
	bool intersectionCircle= false;
	bool intersectionLune = false;
	if(beta <0)
		exit(0);
	else if(beta==0)
		return true;
    	else if(beta <1)
		intersectionCircle = true;
	else
		intersectionLune = true;

	if(this->betaMode == "highDimCircle"){
		if(beta<1)
			beta=1/beta;
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
		
		
             ;
		bool first = true;
		
		std::vector<size_t> neighbors;
		std::vector<std::vector<size_t>> neighborsCircleIntersection;
		for (auto x : simplex){
			double expr1,expr2,expr3;
			std::vector<unsigned> face1;
			face1 = dsimplex;
			face1.erase(std::remove(face1.begin(),face1.end(),x),face1.end());
			std::set<unsigned> face(face1.begin(),face1.end());
			std::vector<double> faceCC ;
			if(face.size()>2)
				faceCC = utils::circumCenter(face,inData.inputData);
			else if(face.size()==2){
				auto first = face.begin();
				std::vector<double> fR;
				std::vector<double> fA = inData.inputData[*first];
      		        	std::advance(first, 1);
				std::vector<double> fB = inData.inputData[*first];
	   			std::transform(fA.begin(), fA.end(), fB.begin(), std::back_inserter(fR),[](double e1,double e2){return ((e1+e2)/2);});
				faceCC = fR;
			}	
			double faceRadius;
   			if(face.size()>2)
   				faceRadius = utils::circumRadius(face,((alphaComplex<alphaNode>*)inData.complex)->distMatrix);
   			else
				faceRadius = pow((*((alphaComplex<alphaNode>*)inData.complex)->distMatrix)[face1[0]][face1[1]]/2,2);
                        
			std::vector<double> hpcoff = utils::nullSpaceOfMatrix(face,inData.inputData,faceCC,sqrt(faceRadius));
			std::vector<double> betaCenter;
			double betaRadius;
             		std::vector<std::vector<double>> betaCenters ;
			bool sameside = false;
			if(intersectionCircle && beta >2){
				if(beta < 3){
					double ratio = sqrt(circumRadius)/sqrt(faceRadius);
		                        betaRadius = sqrt(faceRadius) + (beta-2)*(sqrt(circumRadius)-sqrt(faceRadius));
					betaCenters = utils::betaCentersCalculation(hpcoff, 1+(beta-2)*(ratio-1), sqrt(faceRadius),faceCC);
					
	
				}
				else{

					betaCenters = utils::betaCentersCalculation(hpcoff, beta-1, sqrt(faceRadius),faceCC);
		                 	betaRadius = sqrt(faceRadius)*(beta-1);

                                }
                  		expr1=0;
		  		expr2=0;
                                expr3 =0;
			      	for(unsigned i =0;i<hpcoff.size();i++){
				      	expr1 += hpcoff[i]*circumCenter[i];
					expr2 += hpcoff[i]*betaCenters[0][i];
					expr3 += hpcoff[i]*inData.inputData[x][i];
		      		}
		  		expr1--;
		  		expr2--;
				expr3--;		
				if((expr1> 0 && expr3>0)&&expr2>0||(expr1<0&&expr3<0)&&expr2<0){
					sameside=true;
					betaCenter = betaCenters[1];
					
				}
				else if((expr1>0&&expr3>0)||(expr1<0&&expr3<0)){
					sameside=true;
					if(expr2>0&&expr1>0){
						betaCenter = betaCenters[1];
					}
					else{
						betaCenter = betaCenters[0];
					}
				}else{
					sameside=false;
					if(expr2>0&&expr1>0){
						betaCenter = betaCenters[0];
					}
					else{
						betaCenter = betaCenters[1];
					}
				}
                        }
			else{
				if(face.size()>2){
                  		expr1=0;
		  		expr2=0;
                                expr3 =0;
			      	for(unsigned i =0;i<hpcoff.size();i++){
				      	expr1 += hpcoff[i]*circumCenter[i];
					expr2 += hpcoff[i]*betaCenters[0][i];
					expr3 += hpcoff[i]*inData.inputData[x][i];
		      		}
		  		expr1--;
		  		expr2--;
				expr3--;		
				if((expr1> 0 && expr3>0)&&expr2>0||(expr1<0&&expr3<0)&&expr2<0){
					sameside=true;
				}
				else if((expr1>0&&expr3>0)||(expr1<0&&expr3<0)){
					sameside=true;
				}else{
					sameside=false;
				}
					for(unsigned y =0 ;y< inData.inputData[0].size();y++)
					if(sameside){
						if(intersectionCircle)
							betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);
						else
							betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
					}else{
						if(intersectionCircle)
							betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
						else
							betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);
           				}
				}
				else{
                  		expr1=0;
		  		expr3=0;
			      	for(unsigned i =0;i<hpcoff.size();i++){
				      	expr1 += hpcoff[i]*circumCenter[i];
					expr3 += hpcoff[i]*inData.inputData[x][i];
		      		}
		  		expr1--;
				expr3--;		
				if((expr1>0&&expr3>0)||(expr1<0&&expr3<0)){
					sameside=true;
				}else{
					sameside=false;
				}
				/*	for(unsigned y =0 ;y< inData.inputData[0].size();y++)
					double slope,intercept;
					if(inData.inputData[face1[1]][0]==inData.inputData[face1[0]][0]){
						if(circumCenter[0]<inData.inputData[face1[1]][0]&&inData.inputData[x][0]<inData.inputData[face1[1]][0])
							sameside = true;
	      					else if(circumCenter[0]>inData.inputData[face1[1]][0]&&inData.inputData[x][0]>inData.inputData[face1[1]][0])
		      					sameside = true;
		      				else
	      						sameside=false;
					}
					else{
						slope = (inData.inputData[face1[1]][1] - inData.inputData[face1[0]][1])/(inData.inputData[face1[1]][0]-inData.inputData[face1[0]][0]);
						intercept = (inData.inputData[face1[0]][1]-slope*inData.inputData[face1[0]][0]);
						expr1 = circumCenter[1]-slope*circumCenter[0]-intercept;
			       	  	        expr2 = inData.inputData[x][1]-slope*inData.inputData[x][0]-intercept;
						if(expr1 < 0){
							if(expr2 <0)
								sameside = true;
							else
								sameside = false;
						}
						else{
							if(expr2>0)
								sameside = true;
							else
								sameside = false;
						}
					}
					*/
			
					for(unsigned y =0 ;y< inData.inputData[0].size();y++)
						if(sameside){
							if(intersectionCircle)
								betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);          		
							else
								betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
						}else{
							if(intersectionCircle)
								betaCenter.push_back(beta*circumCenter[y] - (beta-1)*faceCC[y]);
							else
								betaCenter.push_back((2-beta)*circumCenter[y] + (beta-1)*faceCC[y]);
           					}
				}
			}
	   		if(!intersectionCircle || beta<=2){
		     	        betaRadius = utils::vectors_distance(betaCenter,inData.inputData[face1[0]]);
			}
		
		
			std::vector<size_t> neighborsface = tree.neighborhoodIndices(betaCenter, betaRadius); //All neighbors in epsilon-ball
			for(auto x :dsimplex)
		        	neighborsface.erase(std::remove(neighborsface.begin(),neighborsface.end(),x),neighborsface.end());
		       	std::sort (neighborsface.begin(),neighborsface.end());
	      		std::sort (neighbors.begin(),neighbors.end()); 
	    	  	neighborsCircleIntersection.push_back(neighborsface);
	    	
	   		if(!first){
           			 if(intersectionCircle == true){
					std::vector<size_t> v(std::min(neighbors.size(),neighborsface.size()));
					std::vector<size_t>::iterator it;
					it=std::set_intersection (neighbors.begin(), neighbors.end(), neighborsface.begin(), neighborsface.end(), v.begin());
					v.resize(it-v.begin()); 
					neighbors = v;		
				}
				 else{
					std::vector<size_t> v(neighbors.size() + neighborsface.size());
					std::vector<size_t>::iterator it;	
					it=std::set_union (neighbors.begin(), neighbors.end(), neighborsface.begin(), neighborsface.end(), v.begin());
					v.resize(it-v.begin()); 
					neighbors = v;
				}
	    		}
	    		else{
		    		neighbors = neighborsface;
		    		first = false;
	    		}
		}	
		std::vector<size_t> circumneighbors = tree.neighborhoodIndices(circumCenter, sqrt(circumRadius)); //All neighbors in epsilon-ball
		for(auto y :dsimplex){
			circumneighbors.erase(std::remove(circumneighbors.begin(),circumneighbors.end(),y),circumneighbors.end());
		}
		std::sort (circumneighbors.begin(),circumneighbors.end());
		bool circleintersect = false;
		bool brloop = false;
		for(int i=0;i<dsimplex.size()&&!brloop;i++)
			for(int j=i+1;j<dsimplex.size()&&!brloop;j++){			            
				std::vector<size_t> v(std::min(neighborsCircleIntersection[i].size(),neighborsCircleIntersection[j].size()));
				std::vector<size_t>::iterator it;       	
				it=std::set_intersection (neighborsCircleIntersection[i].begin(), neighborsCircleIntersection[i].end(), neighborsCircleIntersection[j].begin(), neighborsCircleIntersection[j].end(), v.begin());
				v.resize(it-v.begin()); 
				if(intersectionCircle)
				{
				std::vector<size_t> newv(std::min(circumneighbors.size(),v.size()));
				std::vector<size_t>::iterator it2;	
				it2=std::set_intersection (circumneighbors.begin(), circumneighbors.end(), v.begin(), v.end(), newv.begin());
				newv.resize(it2-newv.begin());
			        if(newv.size()>0){
					circleintersect = true;
					brloop = true;
				}
				}
				else{
					if(v.size()>0){
						circleintersect = true;
						brloop=true;
					}
				}
			}
                        if(intersectionCircle){
				if(!circleintersect || neighbors.size()==0)
					return true;
				else
					return false;
			}
			else if(!circleintersect)
				return true;
			else
				return false;
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
				if(intersectionLune==true){
					std::vector<size_t> v1(std::min(neighborsfinalLune.size(),neighbors1faces1.size()));
					std::vector<size_t>::iterator it1;
					it1=std::set_intersection (neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1-v1.begin()); 
					neighborsfinalLune = v1;
				}
				else{
					std::vector<size_t> v1(std::max(neighborsfinalLune.size(),neighbors1faces1.size()));
					std::vector<size_t>::iterator it1;				
					it1=std::set_union(neighbors1faces1.begin(), neighbors1faces1.end(), neighborsfinalLune.begin(), neighborsfinalLune.end(), v1.begin());
					v1.resize(it1-v1.begin()); 
					neighborsfinalLune=v1;
				}
			}	
			else
				neighborsfinalLune = neighbors1faces1;
			first= false;
		}
	}
	if(this->betaMode == "highDimLune" && neighborsfinalLune.size() == 0)
		return true;
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
