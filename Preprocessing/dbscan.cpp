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
#include "dbscan.hpp"
#include "utils.hpp"
//////// DBSCAN algorithm for standalone clustering or as an initialization step for DenStream /////////

// basePipe constructor
dbscan::dbscan(){
	procName = "dbscan";
    return;
}
//taking in preprocessor type

std::vector<int> dbscan::cluster(std::vector<std::vector<double>> &data){

   //add separate dbscan path here

}
    // runPipe -> Run the configured functions of this pipeline segment
    pipePacket
    dbscan::runPreprocessor(pipePacket inData)
{ //standalone preprocessor
    // make labels for upscaling
    std::vector<uint_least32_t> upscaleLabels(inData.originalData.size());
    std::iota(std::begin(upscaleLabels), std::end(upscaleLabels), 0);
    ////
    std::cout << "got to outside \n";
    /////////constants//////////
    utils ut;

    //put all data points into the KDTree
    //kdtree -> used for neighbors test

    // while points unprocessed
    //if getNodesinRadius == true (current point, eps, minPts)
    //append to current cluster







          
}


//getNodesinRadius
//query current point in kdtree
// count points near the current points (kd_nearest_range
// if nearest points > minPts 








// configPipe -> configure the function settings of this pipeline segment
bool dbscan::configPreprocessor(std::map<std::string, std::string> configMap){
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
