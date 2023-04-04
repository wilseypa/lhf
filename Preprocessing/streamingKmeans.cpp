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

/**
	@class streamingKmeans
	@brief A class for implementing the streaming k-means algorithm
	@tparam nodeType The type of data points to cluster
*/
// overall goal: get weighted representation of streaming data, then perform k means on that .... Shindler 11

// basePipe constructor

/**
	@brief The overall goal of this code is to get the weighted representation of streaming data and then perform k-means on it.
	@tparam nodeType A template parameter representing the node type of the streaming data.
*/

template<typename nodeType>
/**
	@brief Construct a new streamingKmeans object and set the processor name.
	@tparam nodeType A template parameter representing the node type of the streaming data.
*/
streamingKmeans<nodeType>::streamingKmeans(){
	this->procName = "streamingKmeans";
	return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
/**
	@brief Run the configured functions of this pipeline segment.

	@param nodeType A template parameter representing the node type of the streaming data.

	@param inData A reference to the input data.
*/
template<typename nodeType>
/**
	@brief Run the configured functions of this pipeline segment.
	@param inData A reference to the input data.
*/

/**
	@brief Runs the preprocessor on input data
	@tparam nodeType The datatype of the nodes in the data
	@param inData The input data
 */
void streamingKmeans<nodeType>::runPreprocessor(pipePacket<nodeType> &inData){
	
	if(!this->configured){
		this->ut.writeLog(this->procName,"Preprocessor not configured");
		return;
	}
	
	//Arguments - num_clusters, num_iterations

  int numClusters = 20;

 // int numClusters = 20;

	streamingUtils<nodeType> streamUt;
  //int n = 5;
  int size = inData.workData.size();
	std::cout<< size << "  <-size of data\n";
  int maxFacilities = 20;
//constants used for the online step facility cost increases... increase cost by bounded value as stream progresses (BMORST11 Streaming k-means...)
  const double E = 2.718281828;
  const double alpha = 2.0; //for approx triangle inequality
  const double cofl = 3 * alpha + 2 * (E/(E-1));   //constant online facility location
  const double beta = 2 * alpha * alpha * cofl + (2* alpha);  //constant to increase lower bound
  const double kofl = (6 * alpha) + 1; //constant OFL for upper bound on number of facilities
//doubleFL generates AT MOST kofl(1 + log n)(OPT/L) facilities

 
 std::vector<std::vector<double>> facilities; //(numClusters, std::vector<double>(inData.originalData[0].size(), 0)); //centroids, size = to numClusters "empty set K"
 std::vector<double> omega(std::vector<double>(inData.workData[0].size(), 0)); // vector of values between 0 and 1 based on dim of facilities
 std::vector<double> facilityLabel; // tracks the index of the facilities before they are sorted into the approx facility
 std::vector<double> weight;  //storing weights of points assigned to clusters

for(int d = 0; d<inData.workData[0].size(); d++){  
	omega[d] = randDouble();  // initializing omega (vector that is randomly projected on )
}


//std:: cout<< omega[0] << omega[1] << "  <-omega\n";

std::vector<double> approxFacilities(0);

std::vector<double> approxFacilitiesHat(0);
 std::vector<std::vector<double>> kHat;

/**

	
/**

	@brief Loop through points in a data stream and assign them to clusters, based on delta/f probability.

	@details The method first adds the initial points to the facility set, then loops through the rest of the points

	in the stream, assigning them to the closest cluster. If a point is close enough to none of the existing clusters,

	it is added as a new cluster.

	@param inData The input data to be clustered.

	@param numClusters The number of desired clusters.

	@param maxFacilities The maximum number of facilities that can be added.

	@param omega The random projection matrix.

	@return None.

	@note This method is specific to the kmeans++ algorithm and the data structure used here.
*/

double f = 1/(numClusters*(1 + log(size))); //facility cost f = 1/(k(1+log n))  k clusters, n points, empty set K  //facility==centroid 


std::vector<double> summedClusters(numClusters, 0);

std::vector<double> counts(maxFacilities, 0);

std::vector<int> tempCounts(maxFacilities, 0);
std::vector<unsigned> curLabels;
std::vector< std::pair<double, int>> sortedApproxFacils;
std::vector< std::pair<double, int>> sortedApproxFacilsHat;

std::vector< std::pair<std::vector<double>, int>> clustered;
std::vector< std::pair<std::vector<double>, int>> kHatClustered;


//adding first values to everything so the approximate vectors project correctly

for(unsigned i = 0; i<(numClusters/2); i++){
	facilities.push_back(inData.workData[i]);
	approxFacilities.push_back(dotProd(inData.workData[i], omega));
	sortedApproxFacils.push_back(std::make_pair(approxFacilities[i], i));
	sort(sortedApproxFacils.begin(), sortedApproxFacils.end()); 
	clustered.push_back(std::make_pair(inData.workData[i], i));
	counts[i] ++;
}

int approxSize = sortedApproxFacils.size();

int numFacilities = facilities.size();

int tracker  = 0;

std::vector<std::vector<double>> summedCentroidVectors(numClusters, std::vector<double>(inData.workData[0].size(), 0)); 



//while file stream open
      //while K <= scriptK = klogn (num clusters <= num facilities f) & stream unread
			  //numClusters*log(size)
        
				for(int x = (numClusters/2); x<inData.workData.size()-1; x++) {   //read next point x from the stream
				if (facilities.size() <= maxFacilities){   //first phase goes immediately to reassigning
				 		std::vector<double> y = approxNearestNeighbor(facilities, sortedApproxFacils, omega,  x,  approxSize, pipePacket<nodeType>(inData));
				//		std::cout<< y[0] << "  <-y value \n";
			      double delta = this->ut.vectors_distance(inData.workData[x], y); 	//measure delta = min d(x,y)^2 --> using approx nearest neighbor to get y
			//			std::cout<< inData.originalData[x][0] << inData.originalData[x][1] << " <-current point\n";
			//			std::cout<< delta << " <-delta value \n";
			//			std::cout<< f << " <-f value is \n";
						
				        if(prob(delta/f)){ //as each point arrives either make it a facility or assign it to one based on delta/f prob
                  facilities.push_back(inData.workData[x]);     // K <- K union current point x (add centroid)
				  //        std::cout << "point added to centroids \n";
									approxFacilities.push_back(dotProd(inData.workData[x], omega)); 
									for(int a = 0; a<approxFacilities.size(); a++){//sorting approx facilities ascending so they can be binary searched
										sortedApproxFacils.push_back(std::make_pair(approxFacilities[a], a));
										   
									} 
									sort(sortedApproxFacils.begin(), sortedApproxFacils.end());   // now sortedApproxFacils.first = value, .second = original index
					//				std::cout<< facilities.size() << "  <- facility size \n";
					//				std::cout<< sortedApproxFacils[1].first << sortedApproxFacils[1].second;
								}
				        else{  
					          // add current point to closest cluster center in facilities
						//				std::cout << "point not added to centroids GOT HERE \n";
										double minDist = std::numeric_limits<double>::max();
										for(unsigned c = 0; c<numFacilities; c++){
											// clusterIndex = 0;
											double curDist = this->ut.vectors_distance(inData.workData[x], facilities[c]);
												if(curDist < minDist) {
													minDist = curDist;
													clustered.push_back(std::make_pair(inData.workData[x], c));
													counts[c] ++;
												}
										}
									}
					}  //endif
				
				 else if(facilities.size() > maxFacilities) {   // we reached max facilities count, now evaluate and raise cost
					    f = beta*f; //increasing the cost to add a new centroid 
							for(unsigned q = 0; q<facilities.size(); q++){   //summing the points belonging to each cluster
								if (clustered[q].second == q) {
									for(unsigned dim = 0; dim < facilities[0].size(); dim++){
											summedCentroidVectors[q][dim] += clustered[q].first[dim];  //[q][dim];
									//		std::cout<<clustered[q].first[dim] << "   <-clustered before summed \n";
									//		std::cout<<summedCentroidVectors[q][dim] << " <-summed Centroid vec \n";
									}  
									tempCounts[q] ++; 
							//			std::cout<<summedCentroidVectors[q][0] << " <-summed Centroid vec \n";
								//					std::cout<< "GOT HERE";
								}			
											
							}  
						//	std::cout<< "GOT HERE" ;
							for(int test = 0; test< tempCounts.size()-1; test++){
								if(tempCounts[test] > 0 ){
									tracker++;
								}
					//			std::cout<<tempCounts[test] << " tempCounts \n";
							}    
				//			std::cout<< tempCounts.size() << "  <-tempCount size \n";
			
						//	  std::cout<<summedCentroidVectors.size() << "  <-size of summed centroid vectors \n";
						//		std::cout<<tracker<<"  <-tracker \n";
								summedCentroidVectors.resize(tracker);
						//		  std::cout<<summedCentroidVectors.size() << "  <-size of summed centroid vectors \n";
								tracker = 0; //reset tracker
									for(unsigned i = 0; i < summedCentroidVectors.size(); i++){  //move points x in K to the COM of points for that facility
										for(unsigned dim =0; dim < summedCentroidVectors[0].size(); dim++){
												summedCentroidVectors[i][dim] = summedCentroidVectors[i][dim] / tempCounts[i];
										//		std::cout<<"got here 2"<<" \n";
										
										}
						//				std::cout<<summedCentroidVectors[i][0] << " " << summedCentroidVectors[i][1] << " <-summed before weighting \n"; 
									}  
						
						 weight = counts;   	//set wsubx be number of points assigned to x in K
						 //std::cout<<summedCentroidVectors[0][0] << "  test\n";
						 for(unsigned z = 0; z<summedCentroidVectors.size()-1; z++){
							 		 	approxFacilitiesHat.push_back(dotProd(summedCentroidVectors[z], omega)); 
						 }
						 for(int a = 0; a<approxFacilitiesHat.size(); a++){//sorting approx facilities ascending so they can be binary searched
								sortedApproxFacilsHat.push_back(std::make_pair(approxFacilitiesHat[a], a));
										   
							} 
				
						sort(sortedApproxFacilsHat.begin(), sortedApproxFacilsHat.end()); 
					//		std::cout<<approxFacilitiesHat.size() << "  <-approxHat size \n";
			//			std::cout<<sortedApproxFacilsHat.size() << "  <-sortedHat size \n";
						int sortedHatsize = sortedApproxFacilsHat.size();
						 kHat.push_back(summedCentroidVectors[0]); //initialize K hat containing first facility from K
	              for(int xHat = 0; xHat<summedCentroidVectors.size(); xHat++){  //for each x in K
					//			 		std::cout<< "got to bottom loop" <<std::endl;
										 	 //////// BREAKS RIGHT HERE //////
							       std::vector<double> yHat = approxHat(summedCentroidVectors, sortedApproxFacilsHat, omega,  xHat,   sortedHatsize);
										 //////// BREAKS RIGHT HERE //////
										// 	summedCentroidVectors.clear();
								//			sortedApproxFacilsHat.clear();
							//			 std::cout<< yHat[0] << " <-yhat\n";
							       double deltaHat = this->ut.vectors_distance(summedCentroidVectors[xHat], yHat); 
                         if(prob((weight[xHat]*deltaHat)/f)){   //add old facility to weighted facility set if weight is high enough
									          kHat.push_back(summedCentroidVectors[xHat]);
								 			   }
								 		     else{
									 			   // add x to its closest facility in Khat
															double minDist = std::numeric_limits<double>::max();
														for(unsigned x = 0; x<summedCentroidVectors.size(); x++){
															for(unsigned k = 0; k<kHat.size(); k++){
																	double curDist = this->ut.vectors_distance(summedCentroidVectors[x], kHat[k]);
																		if(curDist < minDist) {
																			minDist = curDist;
																		}
																		for(unsigned d = 0; d<inData.workData[x].size(); d++){
																			kHatClustered.push_back(std::make_pair(summedCentroidVectors[x], k));
																		}
																	
																
																}
														}
							  					}  
					
					//			counts.clear();
								for(unsigned k = 0; k<kHatClustered.size(); k++){   //making Khat = to K 
									facilities.push_back(kHatClustered[k].first);
								}
								kHatClustered.clear();
								for(unsigned r=0; r<facilities.size(); r++ ){    //reassigning approxFacils to Khat
									approxFacilities.push_back(dotProd(facilities[r], omega)); 
										}
								for(int a = 0; a<approxFacilities.size(); a++){
										sortedApproxFacils.push_back(std::make_pair(approxFacilities[a], a));
										   
									} 
									sort(sortedApproxFacils.begin(), sortedApproxFacils.end()); 
								}	//end Khat loop      
						//		std::cout<< x << "  x\n";
											facilities.clear();    //clear everything so main loop iteration can resume
								sortedApproxFacils.clear();
								approxFacilities.clear();
						//		summedCentroidVectors.clear();
								sortedApproxFacilsHat.clear();
								approxFacilitiesHat.clear();
						 
						} //end facility weighting

								
				/*			for(unsigned i = 0; i < summedCentroidVectors.size(); i++){
									//	for(unsigned dim =0; dim < summedCentroidVectors[0].size(); dim++){
											std::cout<<	summedCentroidVectors[i][0] << " ";   // = summedCentroidVectors[i][dim] / counts[i];
											std::cout<<	summedCentroidVectors[i][1]; 
											std::cout << "summed clusters \n";
										}    */
						
						
		
			}  //end reading 	
			for(unsigned c = 0; c<numClusters; c++){
				for(unsigned dim = 0; dim< kHatClustered[c].first.size(); dim++){
						if (std::abs(kHatClustered[c].first[dim] < 4.65825e-20)) {
							kHatClustered[c].first[dim] = 0;
						}

				}
					
			}
				
			std::vector<std::vector<double>> finalClusters;	
			for(unsigned f = 0; f<numClusters; f++){
     				finalClusters.push_back(kHatClustered[f].first);
				}  			

				for(unsigned f = 0; f<finalClusters.size(); f++){
					for(unsigned g = 1; g<finalClusters[0].size(); g++){
							std::cout<<finalClusters[f][g-1];
							std::cout<<"  ";
							std::cout<<finalClusters[f][g];
							 std::cout<<"   <-clustered points \n";
					}
     		 
				}
							
			// now that stream has been exhausted.... Run batch k-means on weighted points K (regular k means)
			//only need to use following 2 kmeans on massive data sets where finalClusters would be  >1000 points... roughly
      //    std::vector<std::vector<double>> clusters =  streamUt.kMeans(finalClusters);	
		//			std::vector<std::vector<double>> ballClusters = streamUt.ballKmeans(clusters); //ball Kmeans to get better centroids




	inData.workData = finalClusters;	
	return;	
}


/**

	@brief Configures the function settings of this pipeline segment.

	@tparam nodeType The data type of the input nodes.

	@param configMap A map containing the configuration parameters and their values.

	@return A boolean value indicating whether the configuration was successful or not.
*/

// configPipe -> configure the function settings of this pipeline segment
template<typename nodeType>
bool streamingKmeans<nodeType>::configPreprocessor(std::map<std::string, std::string> &configMap){
	std::string strDebug;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();
	
	this->ut = utils(strDebug, this->outputFile);
	
	pipe = configMap.find("clusters");
    if(pipe !=configMap.end())
        this->numClusters = std::atoi(configMap["clusters"].c_str());
    else return false;

    pipe = configMap.find("iterations");
	if(pipe != configMap.end())
		this->num_iterations = std::atoi(configMap["iterations"].c_str());
	else return false;	
	
	this->configured = true;
	this->ut.writeDebug(this->procName,"Configured with parameters { clusters: " + configMap["clusters"] + ", iterations: " + configMap["iterations"] + ", debug: " + strDebug + ", outputFile: " + this->outputFile + " }");
	
	return true;
}

/**

	@brief Finds the approximate nearest neighbor for a given data point using the given facilities and projection vector.

	@tparam nodeType The data type of the input nodes.

	@param facilities A vector containing the facilities/centroids.

	@param sortedApproxFacils A vector of pairs containing the inner product of each facility with the projection vector and its index in the facilities vector, sorted in ascending order by the inner product.

	@param omega The projection vector.

	@param x The index of the data point being examined.

	@param size The number of facilities/centroids.

	@param inData The input data packet.

	@return A vector containing the approximate nearest neighbor.
*/

template<typename nodeType>
std::vector<double> streamingKmeans<nodeType>::approxNearestNeighbor(std::vector<std::vector<double>>& facilities, std::vector<std::pair <double, int>>& sortedApproxFacils, std::vector<double> omega, int x, int size, pipePacket<nodeType>(inData)){
  //based on random projection, x is current point being examined, n is number of centroids/facilities
	utils ut;
	double squareDist;
  double projection = dotProd(inData.workData[x], omega);
	int loc = binarySearch(sortedApproxFacils, omega, facilities.size(), projection);
//	std::cout<<loc<<"  <-loc\n";
	//store facilities sorted by their inner product with omega 
	//when new point x arrives, find 2 facilities/centroids that x dotProd omega is between
	//pick which centroid is exactly closer to x , set that equal to delta or "closest facility"
	
	if(loc == -1){
		//nearest centroid is current one
		 squareDist = ut.vectors_distance(inData.workData[x], facilities[sortedApproxFacils[0].second] ); 
//		 std::cout<<squareDist<< " <-squareDist x\n";
		return inData.workData[x];  
	}
	else if (loc == x-1){
		double squareDist =  ut.vectors_distance(inData.workData[x], facilities[sortedApproxFacils[size-1].second]);
//		std::cout<<squareDist<< " <-squareDist x-1\n";
		return facilities[sortedApproxFacils[size-1].second];
	}
 
	
	squareDist = ut.vectors_distance(inData.workData[x], facilities[sortedApproxFacils[loc].second] );
//	std::cout<<squareDist<< " <-squareDist ==0\n";
	double dist = ut.vectors_distance(inData.workData[x], facilities[sortedApproxFacils[loc + 1].second]);
//	std::cout<<dist<< " <-Dist \n";
  if(squareDist <= dist){
		return facilities[sortedApproxFacils[loc].second];
	}
	else {
		squareDist = dist;
		return facilities[sortedApproxFacils[loc + 1].second];
	} 


}

template<typename nodeType>
std::vector<double> streamingKmeans<nodeType>::approxHat(std::vector<std::vector<double>>& summedCentroidVectors, std::vector<std::pair <double, int>>& sortedApproxFacilsHat, std::vector<double> omega, int xHat, int size){
  //based on random projection, x is current point being examined, n is number of centroids/facilities
	utils ut;
	double squareDist;
  double projection = dotProd(summedCentroidVectors[xHat], omega);
	int loc = binarySearch(sortedApproxFacilsHat, omega, size, projection);
//	std::cout<<loc<<"  <-locHat\n";
	//store facilities sorted by their inner product with omega 
	//when new point x arrives, find 2 facilities/centroids that x dotProd omega is between
	//pick which centroid is exactly closer to x , set that equal to delta or "closest facility"
//	std::cout<<"got to approxHat start \n";
	if(loc == -1){
		//nearest centroid is current one
		 squareDist = ut.vectors_distance(summedCentroidVectors[xHat], summedCentroidVectors[sortedApproxFacilsHat[0].second] ); 
//		 std::cout<<squareDist<< " <-squareDist x\n";
//std::cout<<"got to 1\n";
		return summedCentroidVectors[xHat];  
	}
	else if (loc == xHat-1){
		double squareDist =  ut.vectors_distance(summedCentroidVectors[xHat], summedCentroidVectors[sortedApproxFacilsHat[size-1].second]);
//		std::cout<<squareDist<< " <-squareDist x-1\n";
//std::cout<<"got to 2\n";
		return summedCentroidVectors[sortedApproxFacilsHat[size-1].second];
	}
 

	squareDist = ut.vectors_distance(summedCentroidVectors[xHat], summedCentroidVectors[sortedApproxFacilsHat[loc].second] );
//	std::cout<<squareDist<< " <-squareDist ==0\n";

 //std::cout<<sortedApproxFacilsHat[2].first << " sanity check \n";
// std::cout<<summedCentroidVectors[2][0]<< " sanity check \n";
	double dist = ut.vectors_distance(summedCentroidVectors[xHat], summedCentroidVectors[sortedApproxFacilsHat[loc + 1].second]);

//	std::cout<<dist<< " <-Dist \n";
  if(squareDist <= dist){
		return summedCentroidVectors[sortedApproxFacilsHat[loc].second];
	}
	else {
		squareDist = dist;
		return summedCentroidVectors[sortedApproxFacilsHat[loc + 1].second];
	} 

	
}

/**

	This method approximates the closest centroid or facility to a given data point using a random projection-based algorithm.

	@tparam nodeType The data type of each centroid/facility.

	@param summedCentroidVectors A vector of vectors that holds the summed data points for each centroid/facility.

	@param sortedApproxFacilsHat A vector of pairs that holds the inner product of each centroid/facility with a random projection vector and its index in the summedCentroidVectors vector, sorted by the inner product.

	@param omega A vector that holds the random projection vector.

	@param xHat The index of the data point to find the closest centroid/facility to.

	@param size The number of centroids/facilities.

	@return A vector that holds the data point of the closest centroid/facility to the given data point.
*/

template<typename nodeType>
int streamingKmeans<nodeType>::binarySearch(std::vector<std::pair <double, int>>& sorted, std::vector<double> omega, int n, double projection){ //dotProd is target
  //performing binary search on approxFacilities to find starting location for approxNearestNeighbor
//	std::cout <<projection << "projection\n";
//	std::cout << n << " size of approxFacilHat\n";
//	std::cout << sorted[0].first << " target \n";
	if( projection <= sorted[0].first){  //projection = dotProd of current point and omega
		return -1; //if target < dotProd, search unsuccessful
	} 
  if(projection > sorted[n-1].first) {
		return n-1;
	}
	
	int low = 0;
	int high = n-1;
	int mid; 
	//std::cout <<"got here 8\n" ;
	while(high - low >1) {
		mid = round((high+low)/2);
		if(sorted[mid].first>= projection ) {
			high = mid;
		}
		if(sorted[mid].first <= projection) {
			low = mid;
		}		
	}  
return low;


}

!
/**

	@brief Takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1. When new points arrive, dot product is calculated, and 2 centroids x dot omega is between to find nearest facilities faster than calculating nearest neighbor.
	@tparam nodeType
	@param a std::vector of doubles representing the centroid or facility
	@param b std::vector of doubles representing the weight or omega
	@return double representing the dot product
*/

template<typename nodeType>
double streamingKmeans<nodeType>::dotProd(const std::vector<double>& a, const std::vector<double>& b){
	//takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1
	//when new points arrive, dot product calculated, and find 2 centroids x dot omega is between... faster than calc nearest neighbor
	std::vector<double> temp;
 
  std::transform(a.begin(), a.end(), b.begin(),  std::back_inserter(temp), [](double e1, double e2) {return (e1*e2);});
  
  return std::accumulate(temp.begin(), temp.end(), 0.0);
}

/**

	@brief Generates a random double between 0 and 1 using the implementation from stack overflow
	@tparam nodeType
	@return double representing the random number generated
*/

template<typename nodeType>
double streamingKmeans<nodeType>::randDouble(){    // random double between 0 and 1  (stack overflow implementation)
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

/**

	@brief Determines whether an event should happen given a certain probability
	@tparam nodeType
	@param f double representing the probability of an event occurring
	@return true if event should happen, false otherwise
*/

template<typename nodeType>
bool streamingKmeans<nodeType>::prob(double f){
	// returns true with f's probability.
	return streamingKmeans<nodeType>::randDouble() < f;
}

/**

	@brief Generates a random integer between low and high, inclusive.

	@tparam nodeType The data type of the node.

	@param low The lower limit of the range for the random number.

	@param high The upper limit of the range for the random number.

	@return A random integer between low and high, inclusive.

*/

template<typename nodeType>
int streamingKmeans<nodeType>::random(int low, int high){
	return low + ( rand() % (high - low) );
}

/*// configPipe -> configure the function settings of this pipeline segment
bool streamingKmeans::configPreprocessor(std::map<std::string, std::string> configMap){
  auto preprocessor = configMap.find("clusters");
    if(preprocessor !=configMap.end())
        numClusters = std::atoi(configMap["clusters"].c_str());
    else return false;

    preprocessor = configMap.find("iterations");
	if(preprocessor != configMap.end())
		num_iterations = std::atoi(configMap["iterations"].c_str());
	else return false;	
	
	
	return true;
}  */
