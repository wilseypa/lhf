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
#include "dbscan.hpp"
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
utils ut;
dbscan denseUt;


/////////constants//////////
int initPoints = 20; // points to generate p clusters (large data sets, use 1000)
double lambda = 0.25; // decay factor
int epsilon = 16;  //radius of cluster (same as DBSCAN epsilon)
int mu = 10;   //weight of data points in cluster threshold
double beta = 0.2; // outlier threshold
//int pClusterLabel = 0;
int oClusterLabel = 0;
int pClusterIndex;
int oClusterIndex;
double timestamp;
///////initialize p micro clusters... DBSCAN first N points (N has to be less than size of input data to simulate stream)///////
//this returns cluster labels corresponding to current points

std::vector<int>clusterLabels = denseUt.cluster(inData.originalData);


int pClusterLabel = *std::max_element(clusterLabels.begin(), clusterLabels.end()); //start pClusters at highest cluster from DBSCAN
double Tp = (1/lambda)* log((beta*mu)/((beta*mu)-1)); 
//////////// adding points to p or o clusters and updating Tp///// 
double time = 0.0; //placeholder for time measurement
for(int i = initPoints+1; i<inData.originalData.size(); i++){
     time++;
     std::vector<int> tempLabels = merging(inData.originalData, i, clusterLabels, epsilon);  //do merging on point p --> either becomes p micro cluster or o microcluster
     clusterLabels.insert(std::end( clusterLabels), std::begin(tempLabels), std::end(tempLabels));  //append templabels
   if(fmod(time,Tp) == 0){  //prune clusters accordingly
        for(int j = 0; j<clusterLabels.size(); j++){
          for(int c = 1; c<pClusterIndex-1; c++){
            //if weight of current p cluster < beta* mu, delete



          }

        }
        for(int j = 0; j<clusterLabels.size(); j++){
          for(int o = 1; o<oClusterIndex-1; o++){
            double outlierWeightThreshold = (pow(2,((-lambda)*(time -timestamp + Tp))-1))/(pow(2, (-lambda*Tp) -1));  //big epsilon
              //if weight of current o cluster < big epsilon, delete



          }

        }
   }


             
}
  

	return inData;
}


//function merging (takes in point p, p micro clusters, o micro clusters )
std::vector<int> denStream::merging(std::vector<std::vector<double>>& data, int i, std::vector<int> clusterLabels, double epsilon){
    utils ut;
    std::vector<double> dists(data.size(), 0);
    for(int j = 0; j<data.size(); j++){   //checking for nearest p micro cluster to current point
        if(clusterLabels[j] != 0){
           double tempDist = ut.vectors_distance(data[i], data[j]);
           dists[j] = tempDist;
        }
    }
    int minIndex = std::min_element(dists.begin(), dists.end()) - dists.begin();

    if(ut.vectors_distance(data[i], data[minIndex]) < epsilon){   // if (new radius of nearest p microcluster cp < epsilon)
        clusterLabels[i] = 1;       //merge p into cp
    } 
    else{
        for(int j = 0; j<data.size(); j++){   //checking for nearest o micro cluster to current point
            if(clusterLabels[j] == 0){
                double tempDist = ut.vectors_distance(data[i], data[j]);
                dists[j] = tempDist;
            }   
        }

        int minIndex = std::min_element(dists.begin(), dists.end()) - dists.begin();

        if(ut.vectors_distance(data[i], data[minIndex]) < epsilon){   // if (new radius of nearest p microcluster cp < epsilon)
            clusterLabels[i] = 1;       //merge p into co
            //count current cluster index, take weight
            //if weight > beta*mu, remove labels from outlier buffer and make a new pmicro cluster 
        }
        
       //if(new weight of co (w) > Beta*mu)
          //remove co from outlier buffer and create new p micro cluster by co
      //else create a new o micro cluster by p and insert into outlier buffer
    }  

return clusterLabels;
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
