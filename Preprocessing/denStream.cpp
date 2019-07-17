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
int initPoints = 20; // points to generate p clusters (large data sets, use 1000)
double lambda = 0.25; // decay factor
int epsilon = 16;  //radius of cluster (same as DBSCAN epsilon)
int mu = 10;   //weight of data points in cluster threshold
double beta = 0.2; // outlier threshold
int pClusterLabel = 0;
int oClusterLabel = 0;

///////initialize p micro clusters... DBSCAN first N points (N has to be less than size of input data to simulate stream)///////
//this returns cluster labels corresponding to current points



//////////// adding points to p or o clusters and updating Tp///// 
for(int i = initPoints+1; i<inData.workData.originalData.size(); i++){
     double Tp = (1/lambda)* log((beta*mu)/(beta*mu-1)); 
      // get next point
   //do merging on point p --> either becomes p micro cluster or o microcluster


/////////// pruning and cluster maintenance///////

   //if(current time mod Tp = 0)
      //for each p micro cluster
          // if(weight of cp < beta*mu)
                // delete cp
     // for each o micro cluster 
         //update big Epsilon
         //if(weight of co < big Epsilon)
             //delete co
             
}
    
   

  utils ut;

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
