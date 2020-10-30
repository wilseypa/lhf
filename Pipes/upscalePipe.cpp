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



// basePipe constructor
upscalePipe::upscalePipe(){
	pipeType = "Upscale";
	
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
void upscalePipe::runPipe(pipePacket &inData){
	utils ut;
	
	std::vector<std::pair<std::set<unsigned>,std::vector<bettiBoundaryTableEntry>>> upscaleBoundaries;
	
	//Handle two types of boundary sets coming in -
	//		1. Sets with only 2 centroids, indicating a minDist
	//		2. Sets with >2 centroids, indicating a feature for upscaling
	
	std::cout << "Running upscale pipe" << std::endl;
	
	//First, filter and join all of our boundaries by union of sets
	
	for(auto pi = inData.bettiTable.begin(); pi != inData.bettiTable.end(); pi++){
		std::cout << (*pi).bettiDim << ",\t" << (*pi).birth << ",\t" << (*pi).death << ",\t";
		ut.print1DVector((*pi).boundaryPoints);
		
		//Check if this interval lives for longer than scalarV
		if(((*pi).death - (*pi).birth) > scalarV && (*pi).bettiDim > 0){
			
			//Check if there is a set intersection with any of the existing boundary sets
			bool isFound = false;
			int index = 0;
			int firstIntersect = -1;
			
			for(auto bp = upscaleBoundaries.begin(); bp != upscaleBoundaries.end(); bp++){
				if(ut.setIntersect((*bp).first, (*pi).boundaryPoints, true).size() > 0){
					
					if(!isFound){
						upscaleBoundaries[index].first = ut.setUnion((*bp).first, (*pi).boundaryPoints);
						//for(auto bp : pi.boundaryPoints)
						//	upscaleBoundaries[i].insert(bp);
						(*bp).second.push_back((*pi));
						isFound = true;
						firstIntersect = index;
					} else {
						upscaleBoundaries[firstIntersect].first = ut.setUnion((*bp).first, upscaleBoundaries[firstIntersect].first);
						
						for(auto betti : (*bp).second)
							upscaleBoundaries[firstIntersect].second.push_back(betti);
						
						upscaleBoundaries.erase(bp--);
					}
						
				}
				index++;
			}
			if(!isFound){
				std::vector<bettiBoundaryTableEntry> a = {(*pi)};
				upscaleBoundaries.push_back(std::make_pair((*pi).boundaryPoints, a));
			}
			
			inData.bettiTable.erase(pi--);
		}
	}
	
	//Need to re-analyze upscaleBoundaries in case of additional set intersections
	
	std::cout << "Found " << upscaleBoundaries.size() << " independent features to upscale" << std::endl;
	for(auto bound : upscaleBoundaries){
		std::cout <<"\t";
		ut.print1DVector(bound.first);
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
			
			
			auto curwD = pipePacket(subConfigMap,subConfigMap["complexType"]);//args, args["complexType"]);
			
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
				ut.print1DVector(a.boundaryPoints);
				
				//Check if this interval lives for longer than scalarV
				if((a.death - a.birth) > scalarV && a.bettiDim > 0)
					inData.bettiTable.push_back(a);
				
			}
		}
	}
	
	//Merge the new upscaled features into the bettiTable
	
	
	return;
}


void upscalePipe::runSubPipeline(pipePacket& wrData)
{
    if(wrData.workData.size() == 0)
        return;

    outputData(wrData);
    
    

	std::string pipeFuncts = "distMatrix.neighGraph.rips.fast";
    auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;

    //For each '.' separated pipeline function (count of '.' + 1 -> lim)
    for(unsigned i = 0; i < lim; i++)
    {
        auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
        pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

        //Build the pipe component, configure and run
        auto cp = basePipe::newPipe(curFunct, "simplexTree");

        //Check if the pipe was created and configure
        if(cp != 0 && cp->configPipe(subConfigMap))
        {
            //Run the pipe function (wrapper)
            cp->runPipeWrapper(wrData);
        }
        else
        {
            std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
        }
    }


    return;
}




// configPipe -> configure the function settings of this pipeline segment
bool upscalePipe::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;
    subConfigMap = configMap;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = configMap["outputFile"].c_str();
	
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	std::cout << "UPSCALE DIM: " << dim << std::endl;
	
	
	ut = utils(strDebug, outputFile);
	
	//pipe = configMap.find("scalarV");
	//if(pipe != configMap.end())
	//	scalarV = std::atof(configMap["scalarV"].c_str());
	//else return false;
	
	configured = true;
	ut.writeDebug("upscale","Configured with parameters { debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

