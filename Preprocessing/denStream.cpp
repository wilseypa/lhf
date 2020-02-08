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

microCluster::microCluster(int t = 0){
  creationTime = t;
}

void microCluster::insertPoint(std::vector<double> point){

}

std::vector<double> microCluster::getCenter(){
  return CF1; //----------------------------------------------------------------------
}

int microCluster::getCreationTime(){
  return creationTime;
}

double microCluster::getWeight(){
  return weight;
}

double microCluster::getRadius(){

}

// basePipe constructor
denStream::denStream(){
	procName = "DenStream";
    return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
pipePacket denStream::runPreprocessor(pipePacket inData){
  /////////constants//////////
  initPoints = 20; // points to generate p clusters (large data sets, use 1000) //---------------
  minPoints = 5; //For DBSCAN //---------------------------------------------------------------
  lambda = 0.25; // decay factor
  epsilon = 16;  //radius of cluster (same as DBSCAN epsilon)
  mu = 10;   //weight of data points in cluster threshold
  beta = 0.2; // outlier threshold
  streamSpeed = 1; //number of points in one unit time

  Tp = ceil((1/lambda)*log((beta*mu)/((beta*mu)-1)));
  timestamp = 0;
  double xi_den = pow(2, -Tp*lambda)-1; //Denominator for outlier lower weight limit
  int numPerTime = 0;

  ///////initialize p micro clusters... DBSCAN first N points (N has to be less than size of input data to simulate stream)///////
  //this returns cluster labels corresponding to current points
  std::vector<int> clusterLabels = dbscan::cluster(inData.originalData, minPoints, epsilon, initPoints);
  int pClusters = *std::max_element(clusterLabels.begin(), clusterLabels.end()); //Number of clusters found

  std::vector<microCluster> pMicroClusters(pClusters);
  std::vector<microCluster> oMicroClusters;

  //Take labels of each point and insert into appropriate microCluster
  for(int i = 0; i<initPoints; i++){
    if(clusterLabels[i] != -1){ //Not noise
      pMicroClusters[clusterLabels[i]-1].insertPoint(inData.originalData[i]);
    }
  }

  for(int i = initPoints; i<inData.originalData.size(); i++){
    numPerTime++;
    if(numPerTime == streamSpeed){ //1 unit time has elapsed
      numPerTime = 0;
      timestamp++;
    }

    merging(inData.originalData, i);

    if(timestamp % Tp == 0 && numPerTime == 0){
      auto itP = pMicroClusters.begin();
      while(itP != pMicroClusters.end()){
        if(itP->getWeight() < beta*mu){ //Weight too low - no longer a pMicroCluster
          itP = pMicroClusters.erase(itP);
        } else{
          ++itP;
        }
      }

      auto itO = oMicroClusters.begin();
      while(itO != oMicroClusters.end()){
        int T0 = itO->getCreationTime();
        double xi_num = pow(2, -lambda*(timestamp - T0 + Tp)) - 1;
        double xi = xi_num/xi_den; //Outlier lower weight limit

        if(itO->getWeight() < xi){ //Weight too low - no longer an oMicroCluster
          itO = oMicroClusters.erase(itO);
        } else{
          ++itO;
        }
      }
    }

  }

  std::vector<std::vector<double>> centers(pMicroClusters.size());
  for(int i=0; i<pMicroClusters.size(); i++){
    centers[i] = pMicroClusters[i].getCenter();
  }

  inData.originalData = centers;
	return inData;
}

int denStream::nearestPCluster(std::vector<std::vector<double>> &data, int p){
  utils ut;

  int min = 0;
  double minDist = ut.vectors_distance(pMicroClusters[0].getCenter(), data[p]); 

  for(int i=1; i<pMicroClusters.size(); i++){
    double dist = ut.vectors_distance(pMicroClusters[i].getCenter(), data[p]);
    if(dist < minDist){
      minDist = dist;
      min = i;
    }
  }

  return min;
}

int denStream::nearestOCluster(std::vector<std::vector<double>> &data, int p){
  utils ut;

  int min = 0;
  double minDist = ut.vectors_distance(oMicroClusters[0].getCenter(), data[p]); 

  for(int i=1; i<oMicroClusters.size(); i++){
    double dist = ut.vectors_distance(oMicroClusters[i].getCenter(), data[p]);
    if(dist < minDist){
      minDist = dist;
      min = i;
    }
  }

  return min;
}

//function merging (takes in point p, p micro clusters, o micro clusters )
void denStream::merging(std::vector<std::vector<double>> &data, int p){
  int nearestP = nearestPCluster(data, p);
  if(pMicroClusters[nearestP].mergeRadius(data[p]) <= epsilon){
    pMicroClusters[nearestP].insertPoint(data[p]);
    return;
  }

  int nearestO = nearestOCluster(data, p);
  if(oMicroClusters[nearestO].mergeRadius(data[p]) <= epsilon){
    oMicroClusters[nearestO].insertPoint(data[p]);
    if(oMicroClusters[nearestO].getWeight() > beta*mu){
      //--------------------------------- ERASE + add to p
    }
    return;
  }

  microCluster x(timestamp);
  x.insertPoint(data[p]);
  oMicroClusters.push_back(x);
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
