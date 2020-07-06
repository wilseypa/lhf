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
pipePacket upscalePipe::runPipe(pipePacket inData){
	utils ut;
	
	//Handle two types of boundary sets coming in -
	//		1. Sets with only 2 centroids, indicating a minDist
	//		2. Sets with >2 centroids, indicating a feature for upscaling
	
	std::cout << "Running upscale pipe" << std::endl;
	
	for(auto bound : inData.boundaries){
		
		std::cout << bound.size() << std::endl;
		
		if(bound.size() > 2){
			
			
			
			
			
		} else if (bound.size() == 2) {
			//Find minimum distance between 2 indexed partitions
			double minDist = -1;
			std::cout << "Looking for inner cluster distance between " << (*bound.begin()) << " and " << (*bound.rbegin()) << std::endl;
			
			for(unsigned index = 0; index < inData.originalLabels.size(); index++){
				
				std::cout << "a " << inData.originalLabels[index] << std::endl;
				if(inData.originalLabels[index] == (*bound.begin())){
					
					for(unsigned comp = 0; comp < inData.originalLabels.size(); comp++){
						
						if(inData.originalLabels[comp] == (*bound.rbegin())){
							auto dist = ut.vectors_distance(inData.originalData[index], inData.originalData[comp]);
							if(dist == -1 || dist < minDist)
								minDist = dist;			
					
						}
					}
				}
			}
			
			if(minDist > 0)
				std::cout << "Found minDist between two clusters: " << minDist << std::endl;
		
		
		}
		
		
	}
	
	
	
	
	//Merge the new upscaled features into the bettiTable
	
	
	return inData;
}


void upscalePipe::runSubPipeline(pipePacket wrData)
{
    if(wrData.originalData.size() == 0)
        return;

    pipePacket inData = wrData;
    outputData(inData);

	std::string pipeFuncts = "rips.fast";
    auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;

    //For each '.' separated pipeline function (count of '.' + 1 -> lim)
    for(unsigned i = 0; i < lim; i++)
    {
        auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
        pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

        //Build the pipe component, configure and run
        auto *cp = basePipe::newPipe(curFunct, "simplexTree");

        //Check if the pipe was created and configure
        if(cp != 0 && cp->configPipe(subConfigMap))
        {
            //Run the pipe function (wrapper)
            inData = cp->runPipeWrapper(inData);
        }
        else
        {
            std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
        }
    }


    return;
}




// configPipe -> configure the function settings of this pipeline segment
bool upscalePipe::configPipe(std::map<std::string, std::string> configMap){
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
	
	ut = utils(strDebug, outputFile);
	
	pipe = configMap.find("scalarV");
	if(pipe != configMap.end())
		scalarV = std::atof(configMap["scalarV"].c_str());
	else return false;
	
	configured = true;
	ut.writeDebug("upscale","Configured with parameters { debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

