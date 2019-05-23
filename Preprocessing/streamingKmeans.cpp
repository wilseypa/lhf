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
#include "streamingKmeans.hpp"
#include "streamingUtils.hpp"
#include "utils.hpp"
// basePipe constructor
streamingKmeans::streamingKmeans(){
	procName = "streamingKmeans";
	return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket streamingKmeans::runPreprocessor(pipePacket inData){
	//Arguments - num_clusters, num_iterations

	utils ut;
	static std::random_device seed;
	static std::mt19937 gen(seed());
	std::uniform_int_distribution<size_t> distribution(0, inData.workData.originalData.size()-1);

//constants used for the online step (BMORST11 Streaming k-means...)
const float E = 2.718281828;
const float alpha = 2.0; //for approx triangle inequality
const float cofl = 3 * alpha + 2 * (E/(E-1));   //constant online facility location
const float beta = 2 * alpha * alpha * cofl + (2* alpha);  //constant to increase lower bound
const float kofl = (6 * alpha) + 1; //constant OFL for upper bound on number of facilities
// OFL generates AT MOST kofl(1 + log n)(OPT/L) facilities
 std::vector<std::vector<double>> facilities; //centroids
 std::vector<std::vector<double>> omega; // vector of values between 0 and 1 based on dim of facilities









//facility cost f = 1/(k(1+log n))  k clusters, n points, empty set K  //facility==centroid
//as each point arrives either make it a facility or assign it to one based on delta/f prob

//while file stream open
      //while K <= scriptK = klogn (num clusters <= num facilities f) & stream unread
			    //read next point x from the stream


					//measure delta = min d(x,y)^2 --> using approx nearest neighbor

					//if delta/f event occurs
					    // K <- K union current point x
				  // else assign x to closest facility in K

		  // if current stream not exhausted....
			   //while K <= scriptK 
				    //facility <- Beta*facility

						//move points x in K to the COM of points for that facility

						//set wsubx be number of points assigned to x in K

						//intialize K hat containing first facility from K
	              //for each x in K
								   //measure delta = min d(x,y)^2 --> using approx nearest neighbor
									 //if delta/f event occurs 
									     // K hat <- K hat union current point x
									 //else assign x to its closest facility K hat
							  // set K <- K hat
	  // else Run batck k-means on weighted points K
		// Perform ball k-means on set of clusters from previous step
     



	return inData;
}

// configPipe -> configure the function settings of this pipeline segment
bool streamingKmeans::configPreprocessor(std::map<std::string, std::string> configMap){
  auto preprocessor = configMap.find("clusters");
    if(preprocessor !=configMap.end())
        num_clusters = std::atoi(configMap["clusters"].c_str());
    else return false;

    preprocessor = configMap.find("iterations");
	if(preprocessor != configMap.end())
		num_iterations = std::atoi(configMap["iterations"].c_str());
	else return false;	
	
	
	return true;
}

void streamingKmeans:: approxNearestNeighbor(std::vector<std::vector<double>> facilities, float dotProd, int n, float distSquare, int dim ){

}


void streamingKmeans:: binarySearch(std::vector<std::vector<double>> approxFacilities, int n, double target){

}


float dotProd(std::vector<std::vector<double>> facilities, std::vector<std::vector<double>> omega){
	//takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1
	//when new points arrive, dot product calculated, and find 2 centroids x dot omega is between... faster than calc nearest neighbor
	std::vector<double> temp;

  std::transform(facilities.begin(), omega.begin(), std::back_inserter(temp), [](double e1, double e2)  {return e1*e2;});
  
  return std::accumulate(temp.begin(), temp.end(), 0.0);
}