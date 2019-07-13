/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute
 *
 */
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <random>
#include <chrono>
#include <string>
#include <numeric>
#include <iostream>
#include <functional> 
#include <vector>
#include "densityUtils.hpp"
#include "utils.hpp"
//////// DBSCAN algorithm for standalone clustering or as an initialization step for DenStream 

// basePipe constructor
densityUtils::densityUtils(){
	procName = "densityUtils";
    return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket densityUtils::runPreprocessor(pipePacket inData){    //standalone preprocessor
/////////constants//////////
 utils ut;

 //eps - how close points should be to constitute a cluster (need to adjust but preferably small)
 //minPoints - minimum # of points to from a dense region (number of dimensions +1, higher if noisier data/larger data)


	return inData;
}

std::vector<std::vector<double>> densityUtils::dbscan(std::vector<std::vector<double>>& data){   //initialization for DenStream
utils ut; 


return data;

}


















// configPipe -> configure the function settings of this pipeline segment
bool densityUtils::configPreprocessor(std::map<std::string, std::string> configMap){
  /*  auto preprocessor = configMap.find("clusters");
     if(preprocessor !=configMap.end())
        num_clusters = std::atoi(configMap["clusters"].c_str());
    else return false;

    preprocessor = configMap.find("iterations");
	if(preprocessor != configMap.end())
		num_iterations = std::atoi(configMap["iterations"].c_str());
	else return false;  */

	return true;
}
