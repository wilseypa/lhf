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
#include <utility>
#include "streamingKmeans.hpp"
#include "streamingUtils.hpp"
#include "utils.hpp"
// overall goal: get weighted representation of streaming data, then perform k means on that .... Shindler 11

// basePipe constructor
streamingKmeans::streamingKmeans(){
	procName = "streamingKmeans";
	return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket streamingKmeans::runPreprocessor(pipePacket inData){
	//Arguments - num_clusters, num_iterations
  int numClusters = 20;
	utils ut;
	streamingUtils streamUt;
  int n = 5;
  int size = inData.workData.originalData.size();
	std::cout<< size << "  <-size of data\n";
  int maxFacilities = 20;
//constants used for the online step facility cost increases... increase cost by bounded value as stream progresses (BMORST11 Streaming k-means...)
  const double E = 2.718281828;
  const double alpha = 2.0; //for approx triangle inequality
  const double cofl = 3 * alpha + 2 * (E/(E-1));   //constant online facility location
  const double beta = 2 * alpha * alpha * cofl + (2* alpha);  //constant to increase lower bound
  const double kofl = (6 * alpha) + 1; //constant OFL for upper bound on number of facilities
//doubleFL generates AT MOST kofl(1 + log n)(OPT/L) facilities

 
 std::vector<std::vector<double>> facilities; //(numClusters, std::vector<double>(inData.workData.originalData[0].size(), 0)); //centroids, size = to numClusters "empty set K"
 std::vector<double> omega(std::vector<double>(inData.workData.originalData[0].size(), 0)); // vector of values between 0 and 1 based on dim of facilities
 std::vector<double> facilityLabel; // tracks the index of the facilities before they are sorted into the approx facility
 std::vector<double> weight;  //storing weights of points assigned to clusters

for(int d = 0; d<inData.workData.originalData[0].size(); d++){  
	omega[d] = randDouble();  // initializing omega (vector that is randomly projected on )
}
std:: cout<< omega[0] << omega[1] << "  <-omega\n";
std::vector<double> approxFacilities(0);
 std::vector<std::vector<double>> kHat;

double f = 1/(numClusters*(1 + log(size))); //facility cost f = 1/(k(1+log n))  k clusters, n points, empty set K  //facility==centroid 

std::vector<double> summedClusters(num_clusters, 0);
std::vector<std::vector<double>> summedCentroidVectors;
std::vector<double> counts;
std::vector<unsigned> curLabels;
std::vector< std::pair<double, int>> sortedApproxFacils;

//adding first values to everything so the approximate vectors project correctly
facilities.push_back(inData.workData.originalData[0]);
approxFacilities.push_back(dotProd(inData.workData.originalData[0], omega));
sortedApproxFacils.push_back(std::make_pair(approxFacilities[0], 0));
/*std::cout<<sortedApproxFacils[0].first<< "test5\n";
std::cout<<sortedApproxFacils[0].second<< "test4\n";
std::cout<<facilities[sortedApproxFacils[0].second][0]<< "test6 \n";
std::cout<<facilities[0][0] << "test7 \n";  */
int approxSize = sortedApproxFacils.size();
//while file stream open
      //while K <= scriptK = klogn (num clusters <= num facilities f) & stream unread
			  //numClusters*log(size)
        
				for(int x = 1; x<inData.workData.originalData.size(); x++) {   //read next point x from the stream
					int numFacilities = 1;
					while(numFacilities <= maxFacilities){
						unsigned facilityIndex = 1;
				 		std::vector<double> y = approxNearestNeighbor(facilities, sortedApproxFacils, omega,  x,  approxSize, pipePacket(inData));
						std::cout<< y[0] << "  <-y value \n";
			      double delta = ut.vectors_distance(inData.workData.originalData[x], y); 	//measure delta = min d(x,y)^2 --> using approx nearest neighbor to get y
						std::cout<< inData.workData.originalData[x][0] << inData.workData.originalData[x][1] << " <-current point\n";
						std::cout<< delta << " <-delta value \n";
						std::cout<< f << " <-f value is \n";
					

					  
				        if(prob(delta/f)){ //as each point arrives either make it a facility or assign it to one based on delta/f prob
                  facilities.push_back(inData.workData.originalData[x]);     // K <- K union current point x (add centroid)
									facilityLabel.push_back(facilityIndex);
									numFacilities++;
									facilityIndex++;
				          std::cout << "point added to centroids \n";
									approxFacilities.push_back(dotProd(inData.workData.originalData[x], omega)); // need to keep track of sorted approxFacilities
								  for(int a = 1; a<approxFacilities.size(); a++){//sorting approx facilities ascending so they can be binary searched
										sortedApproxFacils.push_back(std::make_pair(approxFacilities[a], a));
										   
									} 
									sort(sortedApproxFacils.begin(), sortedApproxFacils.end());   // now sortedApproxFacils.first = value, .second = original index
									std::cout<< numFacilities << "  <- facility size \n";
									std::cout<< sortedApproxFacils[1].first << sortedApproxFacils[1].second;
								}
						    	
				        else{  
										summedCentroidVectors = facilities;
					          // add current point to closest cluster center in facilities
										std::cout << "point not added to centroids GOT HERE \n";
										unsigned clusterIndex;
										for(unsigned c = 0; c<numFacilities; c++){
											double minDist = std::numeric_limits<double>::max();
											unsigned clusterIndex = 0;
											double curDist = ut.vectors_distance(inData.workData.originalData[x], facilities[c]);
											std::cout<< curDist << "curDist \n";
												if(curDist < minDist) {
													clusterIndex = c;
													minDist = curDist;
												}
										}
												for(unsigned d = 0; d<inData.workData.originalData[x].size(); d++){
														std::cout << "point not added to centroids GOT HERE 3rd loop\n";
														summedCentroidVectors[clusterIndex][d] += inData.workData.originalData[x][d];
												
													}
													curLabels.push_back(clusterIndex);
												counts[clusterIndex] ++;
							//			summedClusters[clusterIndex] += minDist;
									
							  	}
						}
        /*     if(numFacilities = maxFacilities) {   // we reached max facilities count, now evaluate and raise cost
					    f = beta*f; //increasing the cost to add a new centroid
							//move points x in K to the COM of points for that facility
              	for(unsigned i = 0; i < summedCentroidVectors.size(); i++){
										for(unsigned d =0; d < summedCentroidVectors[0].size(); d++){
												summedCentroidVectors[i][d] = summedCentroidVectors[i][d] / counts[i];
										}
								}
						  weight = counts;   	//set wsubx be number of points assigned to x in K
						 kHat.push_back(facilities[0]); //initialize K hat containing first facility from K
	               for(unsigned xHat = 0; xHat<numFacilities; xHat++){  //for each x in K
							       std::vector<double> yHat = approxNearestNeighbor(facilities,  sortedApproxFacils, omega,  x,  size, pipePacket(inData));
							       double deltaHat = ut.vectors_distance(facilities[xHat], yHat); 
                         if(prob(weight[xHat]*deltaHat/f)){
									          kHat.push_back(facilities[xHat]);
								 			   }
								 		     else{
									 			   // add x to its closest facility in Khat
														for(unsigned x = 0; x<facilities.size(); x++){
															for(unsigned k = 0; k<kHat.size(); k++){
																	double minDist = std::numeric_limits<double>::max();
																	unsigned clusterIndex = 0;
																	auto curDist = ut.vectors_distance(facilities[x], kHat[k]);
																		if(curDist < minDist) {
																			clusterIndex = k;
																			minDist = curDist;
																		}
																		for(unsigned d = 0; d<inData.workData.originalData[x].size(); d++){
																			facilities[clusterIndex][d] += inData.workData.originalData[x][d];
																		}
																		curLabels.push_back(clusterIndex);
																		counts[clusterIndex] ++;
																}
														}
							  					}
									}
						 }
						 numFacilities = 0;
						 facilities.clear();
						 facilities = kHat; 
						//		}  //end of phase transition  */
												 
					}   //end of inData loop processing
							 //setting weight adjusted centroids to original centroids		
							
							  // now that stream has been exhausted.... Run batch k-means on weighted points K (regular k means)
      //    std::vector<std::vector<double>> clusters =  streamUt.kMeans(kHat);
		//			std::vector<std::vector<double>> finalClusters = streamUt.ballKmeans(clusters); //ball Kmeans to get better centroids
//inData.workData.originalData = finalClusters;		
	inData.workData.originalData = kHat;	
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

std::vector<double> streamingKmeans::approxNearestNeighbor(std::vector<std::vector<double>> facilities, std::vector< std::pair <double, int> > sortedApproxFacils, std::vector<double> omega, int x, int size, pipePacket(inData)){
  //based on random projection, x is current point being examined, n is number of centroids/facilities
	utils ut;
	double squareDist;
  double projection = dotProd(inData.workData.originalData[x], omega);
	int loc = binarySearch(sortedApproxFacils, omega, facilities.size(), projection);
	std::cout<<loc<<"  <-loc\n";
	//store facilities sorted by their inner product with omega 
	//when new point x arrives, find 2 facilities/centroids that x dotProd omega is between
	//pick which centroid is exactly closer to x , set that equal to delta or "closest facility"
	
	if(loc == -1){
		//nearest centroid is current one
		 squareDist = ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[0].second] ); 
		return inData.workData.originalData[x];  
	}
	else if (loc == x-1){
		double squareDist =  ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[size-1].second]);
		std::cout<<squareDist<< " <-squareDist";
		return facilities[sortedApproxFacils[size-1].second];
	}
 /* else {
			std::cout<<"got here2\n";
			squareDist = ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[loc].second] );
	}  */

	double dist = ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[loc + 1].second]);
  if(squareDist <= dist){
		return facilities[sortedApproxFacils[loc].second];
	}
	else {
		squareDist = dist;
		return facilities[sortedApproxFacils[loc + 1].second];
	} 


}



int streamingKmeans::binarySearch(std::vector< std::pair <double, int> > sortedApproxFacils, std::vector<double> omega, int n, double projection){ //dotProd is target
  //performing binary search on approxFacilities to find starting location for approxNearestNeighbor
	std::cout <<projection << "projection\n";
	std::cout << n << " size of approxFacil\n";
	std::cout << sortedApproxFacils[0].first << " target \n";
	if( projection <= sortedApproxFacils[0].first){  //projection = dotProd of current point and omega
		return -1; //if target < dotProd, search unsuccessful
	} 
  if(projection > sortedApproxFacils[n-1].first) {
		return n-1;
	}
	
	int low = 0;
	int high = n-1;
	int mid; 
	//std::cout <<"got here 8\n" ;
	while(high - low >1) {
		mid = round((high+low)/2);
		if(sortedApproxFacils[mid].first>= projection ) {
			high = mid;
		}
		if(sortedApproxFacils[mid].first <= projection) {
			low = mid;
		}
		else
		{ std::cout<<"approx not big enough";
			return 0;
		}  
		
	}  
return low;


}


//double dotProd(std::vector<std::vector<double>> approxFacilities, std::vector<std::vector<double>> omega){  //d is dimension
double streamingKmeans::dotProd(const std::vector<double>& a, const std::vector<double>& b){
	//takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1
	//when new points arrive, dot product calculated, and find 2 centroids x dot omega is between... faster than calc nearest neighbor
	std::vector<double> temp;
 
  std::transform(a.begin(), a.end(), b.begin(),  std::back_inserter(temp), [](double e1, double e2) {return (e1*e2);});
  
  return std::accumulate(temp.begin(), temp.end(), 0.0);
}

/*double streamingKmeans::dotProd2D(std::vector<std::vector<double>>&  a, std::vector<std::vector<double>> & b){
	//takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1
	//when new points arrive, dot product calculated, and find 2 centroids x dot omega is between... faster than calc nearest neighbor
	std::vector<double> temp;
 
  std::transform(a.begin(), b.begin(), std::back_inserter(temp), [](double e1, double e2)  {return e1*e2;});
  
  return std::accumulate(temp.begin(), temp.end(), 0.0);
}  */

double streamingKmeans::randDouble(){    // random double between 0 and 1  (stack overflow implementation)
	 std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    double currentRandomNumber = unif(rng);
  //      std::cout << currentRandomNumber << std::endl;
    return currentRandomNumber;
} 

bool streamingKmeans::prob(double f){
	// returns true with f's probability.
	return randDouble() < f;
}

int streamingKmeans::random(int low, int high){
	return low + ( rand() % (high - low) );
}



