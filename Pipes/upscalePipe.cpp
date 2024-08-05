/*
 * 
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include "upscalePipe.hpp"
#include "utils.hpp"
#include "pipePacket.hpp"
#include "cluster.hpp"

// basePipe constructor
template <typename nodeType>
upscalePipe<nodeType>::upscalePipe(){
	this->pipeType = "Upscale";
	
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
//
//	Upscale Pipe:
//			1. For each boundary vector identified
//				-Use index to find centroid label
//				-Upscale to original points from label
//				
//		
//
template <typename nodeType>
void upscalePipe<nodeType>::runPipe(pipePacket<nodeType> &inData){

	std::vector<std::pair<std::set<unsigned>, std::vector<bettiBoundaryTableEntry>>> upscaleBoundaries;
	
	//Handle two types of boundary sets coming in -
	//		1. Sets with only 2 centroids, indicating a minDist
	//		2. Sets with >2 centroids, indicating a feature for upscaling
	
	std::cout << "Running upscale pipe" << std::endl;

	//First, filter and join all of our boundaries by union of sets
	
	for(auto pi = inData.bettiTable.begin(); pi != inData.bettiTable.end(); pi++){
		std::cout << pi->bettiDim << ",\t" << pi->birth << ",\t" << pi->death << ",\t";
		this->ut.print1DVector(pi->boundaryPoints);
		
		//Check if this interval lives for longer than scalarV
		if((pi->death - pi->birth) > scalarV && pi->bettiDim > 0){
			
			//Check if there is a set intersection with any of the existing boundary sets
			bool isFound = false;
			auto firstIntersect = upscaleBoundaries.begin();
			
			for(auto bp = upscaleBoundaries.begin(); bp != upscaleBoundaries.end(); bp++){
				if(this->ut.setIntersect(bp->first, pi->boundaryPoints, true).size() > 0){
					
					if(!isFound){
						bp->first = this->ut.setUnion(bp->first, pi->boundaryPoints);
						bp->second.push_back(*pi);

						isFound = true;
						firstIntersect = bp;
					} else {
						firstIntersect->first = this->ut.setUnion(bp->first, firstIntersect->first);
						firstIntersect->second.insert(firstIntersect->second.end(), bp->second.begin(), bp->second.end());
						
						upscaleBoundaries.erase(bp--);
					}
						
				}
			}

			if(!isFound){
				std::vector<bettiBoundaryTableEntry> a = {(*pi)};
				upscaleBoundaries.push_back(std::make_pair(pi->boundaryPoints, a));
			}
			
			inData.bettiTable.erase(pi--);
		}
	}
	
	//Need to re-analyze upscaleBoundaries in case of additional set intersections
	
	std::cout << "Found " << upscaleBoundaries.size() << " independent features to upscale" << std::endl;
	for(auto bound : upscaleBoundaries){
		std::cout <<"\t";
		this->ut.print1DVector(bound.first);
	}
	
	//Upscale each independent boundary
	for(auto bound : upscaleBoundaries){
		
		if(bound.first.size() < 4){
			//Not a suitable boundary for upscaling; revert found features
			for(auto betti : bound.second){
				std::cout << "Reverting betti, boundary too short" << std::endl;
				inData.bettiTable.push_back(betti);
			}
		} else {
			
			auto curwD = pipePacket<nodeType>(subConfigMap, subConfigMap["complexType"]);
			
			for(unsigned index = 0; index < inData.centroidLabels.size(); index++){
				if(bound.first.find(inData.centroidLabels[index]) != bound.first.end()){
					curwD.workData.push_back(inData.inputData[index]);
				}
			}
			
			std::cout << "Gathered " << curwD.workData.size() << " original points" << std::endl;
			
            
            
			runSubPipeline(curwD);
			
			std::cout << std::endl << "_____UPSCALE BETTIS_______" << std::endl;

			for(auto a : curwD.bettiTable){
				std::cout << a.bettiDim << ",\t" << a.birth << ",\t" << a.death << ",\t";
				this->ut.print1DVector(a.boundaryPoints);
				
				//Check if this interval lives for longer than scalarV
				if((a.death - a.birth) > scalarV && a.bettiDim > 0)
					inData.bettiTable.push_back(a);
				
			}
		}
	}
	
	//Merge the new upscaled features into the bettiTable
	
	
	return;
}


template <typename nodeType>
void upscalePipe<nodeType>::runSubPipeline(pipePacket<nodeType>& wrData){
    if(wrData.workData.size() == 0)
        return;

    this->outputData(wrData);
    
    std::cout << "STARTING UPSCALE ROUTINE WITH CLUSTER SIZE: " << wrData.workData.size() << "\t" << this->maxSize << std::endl;

    if(wrData.workData.size() > this->maxSize){
        //Start with the preprocessing function, if enabled
        
        auto clusters = this->maxSize;
        
        auto algo = kmeansplusplus();
        algo.clusterData(wrData.workData, wrData.workData, wrData.centroidLabels, clusters, 100, -1);
        
        
        auto sv = subConfigMap.find("scalarV");
        if(sv == subConfigMap.end()){
            auto scalar = std::atof(subConfigMap["scalar"].c_str());
            subConfigMap["scalarV"] = std::to_string(scalar * utils::computeMaxRadius(clusters, wrData.workData, wrData.inputData, wrData.centroidLabels));
            std::cout << "Using scalarV: " << subConfigMap["scalarV"] << std::endl;
        }
        
        
        
    }


	std::string pipeFuncts = "distMatrix.neighGraph.incrementalPersistence";
    auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;

    //For each '.' separated pipeline function (count of '.' + 1 -> lim)
    for(unsigned i = 0; i < lim; i++){
        auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
        pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

        //Build the pipe component, configure and run
        auto cp = basePipe<nodeType>::newPipe(curFunct, subConfigMap["complexType"]);

        //Check if the pipe was created and configure
        if(cp != 0 && cp->configPipe(subConfigMap)){
            //Run the pipe function (wrapper)
            cp->runPipeWrapper(wrData);
        } else{
            std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
        }
    }


    return;
}




// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool upscalePipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;
    subConfigMap = configMap;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();
        
	pipe = configMap.find("maxSize");
	if(pipe != configMap.end())
		this->maxSize = std::atoi(configMap["maxSize"].c_str());
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	std::cout << "UPSCALE DIM: " << dim << std::endl;
	
	
	this->ut = utils(strDebug, this->outputFile);
	
	pipe = configMap.find("scalarV");
	if(pipe != configMap.end())
		this->scalarV = std::atof(configMap["scalarV"].c_str());
	else{
		this->ut.writeDebug("upscale","No scalarV set; this indicates no partitioning was done and upscaling is not required");
		return false;
	}
	this->configured = true;
	this->ut.writeDebug("upscale","Configured with parameters { debug: " + strDebug + ", outputFile: " + this->outputFile + " }");
	
	return true;
}

template class upscalePipe<simplexNode>;
template class upscalePipe<alphaNode>;
template class upscalePipe<witnessNode>;
