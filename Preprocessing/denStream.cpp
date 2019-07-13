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
#include "denStream.hpp"
#include "utils.hpp"
/////// Based off algorithm outlined in Cao et al 06 "Density-Based Clustering over an Evolving Data Stream with Noise"/////

// basePipe constructor
denStream::denStream(){
	procName = "DenStream";
    return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket denStream::runPreprocessor(pipePacket inData){
/////////constants//////////





///////initialize p micro clusters... DBSCAN first N points??///////

//////////// adding points to p or o clusters and updating Tp///// 
    //set Tp = (1/lambda)* log((beta*mu)/(beta*mu-1))
    // get next point
   //do merging on point p --> either becomes p mc or o mc








/////////// pruning and cluster maintenance///////

   //if(current time mod Tp = 0)
      //for each p micro cluster
          // if(weight of cp < beta*mu)
                // delete cp
     // for each o micro cluster 
         //update big Epsilon
         //if(weight of co < big Epsilon)
             //delete co
             
  










  utils ut;

	return inData;
}






// configPipe -> configure the function settings of this pipeline segment
bool denStream::configPreprocessor(std::map<std::string, std::string> configMap){
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
