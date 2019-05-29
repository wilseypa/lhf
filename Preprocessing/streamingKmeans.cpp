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
  int numClusters;
	utils ut;
	static std::random_device seed;
	static std::mt19937 gen(seed());
	std::uniform_int_distribution<size_t> distribution(0, inData.workData.originalData.size()-1);
int n = 5;
int size = inData.workData.originalData.size();
//constants used for the online step (BMORST11 Streaming k-means...)
const float E = 2.718281828;
const float alpha = 2.0; //for approx triangle inequality
const float cofl = 3 * alpha + 2 * (E/(E-1));   //constant online facility location
const float beta = 2 * alpha * alpha * cofl + (2* alpha);  //constant to increase lower bound
const float kofl = (6 * alpha) + 1; //constant OFL for upper bound on number of facilities

// OFL generates AT MOST kofl(1 + log n)(OPT/L) facilities
 std::vector<std::vector<double>> facilities; //centroids, size = to numClusters "empty set K"
 std::vector<double> omega; // vector of values between 0 and 1 based on dim of facilities
 std::vector<double> facilityLabel; // tracks the index of the facilities before they are sorted into the approx facility

for(int d = 0; d<inData.workData.originalData[0].size(); d++){
	omega[d] = randDouble();  // initializing omega

}
std::vector<std::vector<double>> approxFacilities;

//facility cost f = 1/(k(1+log n))  k clusters, n points, empty set K  //facility==centroid 
float f = 1/(numClusters*(1+ log(inData.workData.originalData.size())));
//as each point arrives either make it a facility or assign it to one based on delta/f prob

//while file stream open
      //while K <= scriptK = klogn (num clusters <= num facilities f) & stream unread
			while(facilities.size() <= numClusters*log(inData.workData.originalData.size()) ) {
				for(int x = 0; x<inData.workData.originalData.size(); x++) {   //read next point x from the stream
			//	 double delta = approxNearestNeighbor(inData.workData.originalData[x], omega[x], );
				
					//measure delta = min d(x,y)^2 --> using approx nearest neighbor

					//if delta/f event occurs
					    // K <- K union current point x
				  // else assign x to closest facility in K
				
				
				
				
				
				
				
				}
			}
			   


				

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

std::vector<double> streamingKmeans:: approxNearestNeighbor(std::vector<std::vector<double>> facilities, std::vector<double> approxFacilities, std::vector<double> facilityLabel, std::vector<double> omega, int x, int size, pipePacket(inData)){
  //based on random projection 
	//x is current point being examined
	//n is number of centroids/facilities
	utils ut;
  double projection = dotProd(inData.workData.originalData[x], omega);
	int loc = binarySearch(approxFacilities, omega, size, projection);
	//store facilities sorted by their inner product with omega 
	//when new point x arrives, find 2 facilities/centroids that x dotProd omega is between
	//pick which centroid is exactly closer to x , set that equal to delta or "closest facility"
	if(loc == -1){
		//nearest centroid is current
		double squareDist = ut.vectors_distance(inData.workData.originalData[x], facilities[facilityLabel[0]] ); //dist between current point and approxFacil corresponding to current centroid denoted by facilTracker
		
	}
	else if (loc == x-1){
		//nearest is approx[n-1]
		double squareDist =  ut.vectors_distance(inData.workData.originalData[x], facilities[facilityLabel[size-1]]);
	}

	double squareDist = ut.vectors_distance(inData.workData.originalData[x], facilities[facilityLabel[loc]] );

	double dist = ut.vectors_distance(inData.workData.originalData[x], facilities[facilityLabel[loc + 1]]);
  if(squareDist <= dist){
		return facilities[facilityLabel[loc]];
	}
	else {
		squareDist = dist;
		return facilities[facilityLabel[loc + 1]];
	}











	return;

}


int streamingKmeans:: binarySearch(std::vector<double> approxFacilities, std::vector<double> omega, int n, double target){ //dotProd is target
  //performing binary search on approxFacilities to find starting location for approxNearestNeighbor
	if( target < approxFacilities[0]){  //projection = dotProd of 
		return -1; //if target < dotProd, search unsuccessful
	} 
  if(target > approxFacilities[n-1]) {
		return n-1;
	}
	int low = 0;
	int high = n-1;
	int mid; 
	while(high - low >1) {
		mid = (high+low)/2;
		if(approxFacilities[mid]>= target ) {
			high = mid;
		}
		if(approxFacilities[mid] <= target) {
			low = mid;
		}
	}
return low;


}


//double dotProd(std::vector<std::vector<double>> approxFacilities, std::vector<std::vector<double>> omega){  //d is dimension
double streamingKmeans::dotProd(const std::vector<double>& a, const std::vector<double>& b){
	//takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1
	//when new points arrive, dot product calculated, and find 2 centroids x dot omega is between... faster than calc nearest neighbor
	std::vector<double> temp;
 
  std::transform(a.begin(), b.begin(), std::back_inserter(temp), [](double e1, double e2)  {return e1*e2;});
  
  return std::accumulate(temp.begin(), temp.end(), 0.0);
}

double streamingKmeans::dotProd2D(std::vector<std::vector<double>>&  a, std::vector<std::vector<double>> & b){
	//takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1
	//when new points arrive, dot product calculated, and find 2 centroids x dot omega is between... faster than calc nearest neighbor
	std::vector<double> temp;
 
  std::transform(a.begin(), b.begin(), std::back_inserter(temp), [](double e1, double e2)  {return e1*e2;});
  
  return std::accumulate(temp.begin(), temp.end(), 0.0);
}

double streamingKmeans::randDouble(){
    return rand()/(double(RAND_MAX)+1);
} 

bool streamingKmeans::prob(double f){
	// returns true with f's probability.
	return randDouble() < f;
}

int streamingKmeans::random(int low, int high){
	return low + ( rand() % (high - low) );
}
