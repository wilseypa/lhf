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
  //int n = 5;
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
std::vector<double> approxFacilitiesHat(0);
 std::vector<std::vector<double>> kHat;

double f = 1/(numClusters*(1 + log(size))); //facility cost f = 1/(k(1+log n))  k clusters, n points, empty set K  //facility==centroid 

std::vector<double> summedClusters(num_clusters, 0);

std::vector<double> counts(maxFacilities, 0);
std::vector<int> tempCounts(maxFacilities, 0);
std::vector<unsigned> curLabels;
std::vector< std::pair<double, int>> sortedApproxFacils;
std::vector< std::pair<double, int>> sortedApproxFacilsHat;
std::vector< std::pair<std::vector<double>, int>> clustered;
std::vector< std::pair<std::vector<double>, int>> kHatClustered;


//adding first values to everything so the approximate vectors project correctly

for(unsigned i = 0; i<(numClusters/2); i++){
	facilities.push_back(inData.workData.originalData[i]);
	approxFacilities.push_back(dotProd(inData.workData.originalData[i], omega));
	sortedApproxFacils.push_back(std::make_pair(approxFacilities[i], i));
	sort(sortedApproxFacils.begin(), sortedApproxFacils.end()); 
	clustered.push_back(std::make_pair(inData.workData.originalData[i], i));
	counts[i] ++;
}
int approxSize = sortedApproxFacils.size();
int numFacilities = facilities.size();
int tracker  = 0;
std::vector<std::vector<double>> summedCentroidVectors(numClusters, std::vector<double>(inData.workData.originalData[0].size(), 0)); 

/*std::cout<<inData.workData.originalData[x][0];
std::cout<< " test point\n" <<std::endl;
std::cout<<inData.workData.originalData[x][1];
std::cout<< " test point\n" <<std::endl;
std::vector<double> y = approxNearestNeighbor(facilities, sortedApproxFacils, omega,  x,  approxSize, pipePacket(inData));
std::cout << y[0] << "   <-y ]\n";  */
   

//while file stream open
      //while K <= scriptK = klogn (num clusters <= num facilities f) & stream unread
			  //numClusters*log(size)
        
				for(int x = (numClusters/2); x<inData.workData.originalData.size(); x++) {   //read next point x from the stream
				if (facilities.size() <= maxFacilities){   //first phase goes immediately to reassigning
				 		std::vector<double> y = approxNearestNeighbor(facilities, sortedApproxFacils, omega,  x,  approxSize, pipePacket(inData));
						std::cout<< y[0] << "  <-y value \n";
			      double delta = ut.vectors_distance(inData.workData.originalData[x], y); 	//measure delta = min d(x,y)^2 --> using approx nearest neighbor to get y
						std::cout<< inData.workData.originalData[x][0] << inData.workData.originalData[x][1] << " <-current point\n";
						std::cout<< delta << " <-delta value \n";
						std::cout<< f << " <-f value is \n";
						
				        if(prob(delta/f)){ //as each point arrives either make it a facility or assign it to one based on delta/f prob
                  facilities.push_back(inData.workData.originalData[x]);     // K <- K union current point x (add centroid)
								//	numFacilities++;
				          std::cout << "point added to centroids \n";
						//			double dotProductTest = dotProd(inData.workData.originalData[x], omega);
									approxFacilities.push_back(dotProd(inData.workData.originalData[x], omega)); 
									for(int a = 0; a<approxFacilities.size(); a++){//sorting approx facilities ascending so they can be binary searched
										sortedApproxFacils.push_back(std::make_pair(approxFacilities[a], a));
										   
									} 
									sort(sortedApproxFacils.begin(), sortedApproxFacils.end());   // now sortedApproxFacils.first = value, .second = original index
									std::cout<< facilities.size() << "  <- facility size \n";
									std::cout<< sortedApproxFacils[1].first << sortedApproxFacils[1].second;
								}
				        else{  
					          // add current point to closest cluster center in facilities
										std::cout << "point not added to centroids GOT HERE \n";
										double minDist = std::numeric_limits<double>::max();
										for(unsigned c = 0; c<numFacilities; c++){
											// clusterIndex = 0;
											double curDist = ut.vectors_distance(inData.workData.originalData[x], facilities[c]);
												if(curDist < minDist) {
													minDist = curDist;
													clustered.push_back(std::make_pair(inData.workData.originalData[x], c));
													counts[c] ++;
												}
										}
									}
					}  //endif
				
				 else if(facilities.size() > maxFacilities) {   // we reached max facilities count, now evaluate and raise cost
					    f = beta*f; //increasing the cost to add a new centroid
					//	  std::vector<std::vector<double>> summedCentroidVectors(numClusters, std::vector<double>(inData.workData.originalData[0].size(), 0)); 
							for(unsigned q = 0; q<facilities.size(); q++){   //summing the points belonging to each cluster
								if (clustered[q].second == q) {
									for(unsigned dim = 0; dim < facilities[0].size(); dim++){
											summedCentroidVectors[q][dim] += clustered[q].first[dim];  //[q][dim];
									//		std::cout<<clustered[q].first[dim] << "   <-clustered before summed \n";
									//		std::cout<<summedCentroidVectors[q][dim] << " <-summed Centroid vec \n";
									}  
									tempCounts[q] ++; 
										std::cout<<summedCentroidVectors[q][0] << " <-summed Centroid vec \n";
													std::cout<< "GOT HERE";
								}			
											
							}  
						//	std::cout<< "GOT HERE" ;
							for(int test = 0; test< tempCounts.size()-1; test++){
								if(tempCounts[test] > 0 ){
									tracker++;
								}
								std::cout<<tempCounts[test] << " tempCounts \n";
							}    
							std::cout<< tempCounts.size() << "  <-tempCount size \n";
					//		tracker = 6;
							  std::cout<<summedCentroidVectors.size() << "  <-size of summed centroid vectors \n";
								std::cout<<tracker<<"  <-tracker \n";
								summedCentroidVectors.resize(tracker);
								  std::cout<<summedCentroidVectors.size() << "  <-size of summed centroid vectors \n";
								tracker = 0;
									for(unsigned i = 0; i < summedCentroidVectors.size(); i++){  //move points x in K to the COM of points for that facility
										for(unsigned dim =0; dim < summedCentroidVectors[0].size(); dim++){
												summedCentroidVectors[i][dim] = summedCentroidVectors[i][dim] / tempCounts[i];
										//		std::cout<<"got here 2"<<" \n";
										
										}
										std::cout<<summedCentroidVectors[i][0] << " " << summedCentroidVectors[i][1] << " <-summed before weighting \n"; 
									}  
						
						 weight = counts;   	//set wsubx be number of points assigned to x in K
						 std::cout<<summedCentroidVectors[0][0] << "  test\n";
						 for(unsigned z = 0; z<summedCentroidVectors.size()-1; z++){
							 		 	approxFacilitiesHat.push_back(dotProd(summedCentroidVectors[z], omega)); 
						 }
						 for(int a = 0; a<approxFacilitiesHat.size(); a++){//sorting approx facilities ascending so they can be binary searched
								sortedApproxFacilsHat.push_back(std::make_pair(approxFacilitiesHat[a], a));
										   
							} 
				
						sort(sortedApproxFacilsHat.begin(), sortedApproxFacilsHat.end()); 
							std::cout<<approxFacilitiesHat.size() << "  <-approxHat size \n";
						std::cout<<sortedApproxFacilsHat.size() << "  <-sortedHat size \n";
						 kHat.push_back(summedCentroidVectors[0]); //initialize K hat containing first facility from K
	              for(int xHat = 0; xHat<summedCentroidVectors.size(); xHat++){  //for each x in K
								 		std::cout<< "got to bottom loop" <<std::endl;
										 	 //////// BREAKS RIGHT HERE //////
							       std::vector<double> yHat = approxHat(kHat, sortedApproxFacilsHat, omega,  xHat,  sortedApproxFacilsHat.size());
										 //////// BREAKS RIGHT HERE //////
										 	summedCentroidVectors.clear();
											sortedApproxFacilsHat.clear();
										 std::cout<< yHat[0] << " <-yhat\n";
							       double deltaHat = ut.vectors_distance(summedCentroidVectors[xHat], yHat); 
                         if(prob((weight[xHat]*deltaHat)/f)){   //add old facility to weighted facility set if weight is high enough
									          kHat.push_back(summedCentroidVectors[xHat]);
								 			   }
								 		     else{
									 			   // add x to its closest facility in Khat
															double minDist = std::numeric_limits<double>::max();
														for(unsigned x = 0; x<summedCentroidVectors.size(); x++){
															for(unsigned k = 0; k<kHat.size(); k++){
																	double curDist = ut.vectors_distance(summedCentroidVectors[x], kHat[k]);
																		if(curDist < minDist) {
																			minDist = curDist;
																		}
																		for(unsigned d = 0; d<inData.workData.originalData[x].size(); d++){
																			//facilities[k][d] += inData.workData.originalData[x][d];
																			kHatClustered.push_back(std::make_pair(summedCentroidVectors[x], k));
																		}
																	
																	//	counts[k] ++;
																}
														}
							  					}  
								facilities.clear();    //clear everything so main loop iteration can resume
								sortedApproxFacils.clear();
								approxFacilities.clear();
								summedCentroidVectors.clear();
								sortedApproxFacilsHat.clear();
								approxFacilitiesHat.clear();
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
						 
						} //end facility weighting

								
				/*			for(unsigned i = 0; i < summedCentroidVectors.size(); i++){
									//	for(unsigned dim =0; dim < summedCentroidVectors[0].size(); dim++){
											std::cout<<	summedCentroidVectors[i][0] << " ";   // = summedCentroidVectors[i][dim] / counts[i];
											std::cout<<	summedCentroidVectors[i][1]; 
											std::cout << "summed clusters \n";
										}    */
						
						
		
			}  //end reading 					
							
						
           
				//		 numFacilities = 0;
				//		 facilities.clear();

						//		}  //end of phase transition  */
												 
					  // } end of inData loop processing
							 //setting weight adjusted centroids to original centroids		
							
							  // now that stream has been exhausted.... Run batch k-means on weighted points K (regular k means)
      //    std::vector<std::vector<double>> clusters =  streamUt.kMeans(kHat);
		//			std::vector<std::vector<double>> finalClusters = streamUt.ballKmeans(clusters); //ball Kmeans to get better centroids
//inData.workData.originalData = finalClusters;	
/*  std::vector<std::vector<double>> finalClusters;
	for(unsigned f = 0; f<numClusters; f++){
     finalClusters.push_back(kHatClustered[f].first);

	}  */

	/*		for(unsigned f = 0; f<finalClusters.size(); f++){
				for(unsigned g = 0; g<finalClusters[0].size(); g++){
						std::cout<<finalClusters[f][g];
				}
     				std::cout<<"   <-clustered points\n";

	}  */

/*	for(unsigned f = 0; f<clustered.size(); f++){
				for(unsigned g = 0; g<facilities[0].size(); g++){
						std::cout<<clustered[f].first[g];
					
				}
					std::cout<<"   <-clustered points";
					std::cout<<clustered[f].second <<  "   <- idx \n";
     			

	}  */
 	  



 // kHatClustered.first = facilities;	
	inData.workData.originalData = facilities;	
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

std::vector<double> streamingKmeans::approxNearestNeighbor(std::vector<std::vector<double>>& facilities, std::vector<std::pair <double, int>>& sortedApproxFacils, std::vector<double> omega, int x, int size, pipePacket(inData)){
  //based on random projection, x is current point being examined, n is number of centroids/facilities
	utils ut;
	double squareDist;
  double projection = dotProd(inData.workData.originalData[x], omega);
	int loc = binarySearch(sortedApproxFacils, omega, facilities.size(), projection);
//	std::cout<<loc<<"  <-loc\n";
	//store facilities sorted by their inner product with omega 
	//when new point x arrives, find 2 facilities/centroids that x dotProd omega is between
	//pick which centroid is exactly closer to x , set that equal to delta or "closest facility"
	
	if(loc == -1){
		//nearest centroid is current one
		 squareDist = ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[0].second] ); 
//		 std::cout<<squareDist<< " <-squareDist x\n";
		return inData.workData.originalData[x];  
	}
	else if (loc == x-1){
		double squareDist =  ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[size-1].second]);
//		std::cout<<squareDist<< " <-squareDist x-1\n";
		return facilities[sortedApproxFacils[size-1].second];
	}
 
	
	squareDist = ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[loc].second] );
//	std::cout<<squareDist<< " <-squareDist ==0\n";
	double dist = ut.vectors_distance(inData.workData.originalData[x], facilities[sortedApproxFacils[loc + 1].second]);
//	std::cout<<dist<< " <-Dist \n";
  if(squareDist <= dist){
		return facilities[sortedApproxFacils[loc].second];
	}
	else {
		squareDist = dist;
		return facilities[sortedApproxFacils[loc + 1].second];
	} 


}

std::vector<double> streamingKmeans::approxHat(std::vector<std::vector<double>>& kHat, std::vector<std::pair <double, int>>& sortedApproxFacils, std::vector<double> omega, int xHat, int size){
  //based on random projection, x is current point being examined, n is number of centroids/facilities
	utils ut;
	double squareDist;
  double projection = dotProd(kHat[xHat], omega);
	int loc = binarySearch(sortedApproxFacils, omega, kHat.size(), projection);
//	std::cout<<loc<<"  <-loc\n";
	//store facilities sorted by their inner product with omega 
	//when new point x arrives, find 2 facilities/centroids that x dotProd omega is between
	//pick which centroid is exactly closer to x , set that equal to delta or "closest facility"
	
	if(loc == -1){
		//nearest centroid is current one
		 squareDist = ut.vectors_distance(kHat[xHat], kHat[sortedApproxFacils[0].second] ); 
//		 std::cout<<squareDist<< " <-squareDist x\n";
		return kHat[xHat];  
	}
	else if (loc == xHat-1){
		double squareDist =  ut.vectors_distance(kHat[xHat], kHat[sortedApproxFacils[size-1].second]);
//		std::cout<<squareDist<< " <-squareDist x-1\n";
		return kHat[sortedApproxFacils[size-1].second];
	}
 
	
	squareDist = ut.vectors_distance(kHat[xHat], kHat[sortedApproxFacils[loc].second] );
//	std::cout<<squareDist<< " <-squareDist ==0\n";
	double dist = ut.vectors_distance(kHat[xHat], kHat[sortedApproxFacils[loc + 1].second]);
//	std::cout<<dist<< " <-Dist \n";
  if(squareDist <= dist){
		return kHat[sortedApproxFacils[loc].second];
	}
	else {
		squareDist = dist;
		return kHat[sortedApproxFacils[loc + 1].second];
	} 


}



int streamingKmeans::binarySearch(std::vector<std::pair <double, int>>& sortedApproxFacils, std::vector<double> omega, int n, double projection){ //dotProd is target
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
	}  
return low;


}



double streamingKmeans::dotProd(const std::vector<double>& a, const std::vector<double>& b){
	//takes dot product of facilities centroids and omega, where omega is d dimensions large uniformly distributed between 0,1
	//when new points arrive, dot product calculated, and find 2 centroids x dot omega is between... faster than calc nearest neighbor
	std::vector<double> temp;
 
  std::transform(a.begin(), a.end(), b.begin(),  std::back_inserter(temp), [](double e1, double e2) {return (e1*e2);});
  
  return std::accumulate(temp.begin(), temp.end(), 0.0);
}

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


