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
	
	std::vector<std::set<unsigned>> upscaleBoundaries;
	
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
			
			for(unsigned i = 0; i < upscaleBoundaries.size(); i++){
				if(ut.setIntersect(upscaleBoundaries[i], (*pi).boundaryPoints, true).size() > 0){
					upscaleBoundaries[i] = ut.setUnion(upscaleBoundaries[i], (*pi).boundaryPoints);
					//for(auto bp : pi.boundaryPoints)
					//	upscaleBoundaries[i].insert(bp);
					isFound = true;
					break;
				}
				
			}
			if(!isFound){
				upscaleBoundaries.push_back((*pi).boundaryPoints);
			}
			
			inData.bettiTable.erase(pi--);
		}
	}
	
	//Need to re-analyze upscaleBoundaries in case of additional set intersections
	
	//Upscale each independent boundary
	for(auto bound : upscaleBoundaries){
		
		
		auto curwD = pipePacket(subConfigMap,subConfigMap["complexType"]);//args, args["complexType"]);
		
		for(unsigned index = 0; index < inData.originalLabels.size(); index++){
			
			if(bound.find(inData.originalLabels[index]) != bound.end())
				curwD.originalData.push_back(inData.originalData[index]);
		
		}
		
		std::cout << "Gathered " << curwD.originalData.size() << " original points" << std::endl;
		
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
	
	//Merge the new upscaled features into the bettiTable
	
	
	return;
}


void upscalePipe::runSubPipeline(pipePacket& wrData)
{
    if(wrData.originalData.size() == 0)
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
	
	ut = utils(strDebug, outputFile);
	
	//pipe = configMap.find("scalarV");
	//if(pipe != configMap.end())
	//	scalarV = std::atof(configMap["scalarV"].c_str());
	//else return false;
	
	configured = true;
	ut.writeDebug("upscale","Configured with parameters { debug: " + strDebug + ", outputFile: " + outputFile + " }");
	
	return true;
}

