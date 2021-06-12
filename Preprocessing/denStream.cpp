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

microCluster::microCluster(int t, int dim, double l){
  creationTime = t;
  lambda = l;
  center = std::vector<double>(dim, 0); //Maintain center (CF1/w) instead of sum of points (CF1) to prevent overflow
  avgSquares = std::vector<double>(dim, 0); //Maintain average of squares (CF2/w) instead sum of squares (CF2) to prevent overflow
  weight = 0;
  editTime = 0; //Maintain last edited time to decay weight correctly
}

void microCluster::insertPoint(std::vector<double> point, int timestamp){
  double newWeight = weight*pow(2, -lambda*(timestamp - editTime)) + 1; //Weight decays before new point is inserted

  if(weight == 0){ //No points in cluster
    center = point;
    for(int i=0; i<point.size(); i++){
      avgSquares[i] = point[i]*point[i];
    }    
  } else{ //Add point to cluster - update center and average sum of squares
    for(int i=0; i<point.size(); i++){
      center[i] = center[i] + (point[i]-center[i])/newWeight;
      avgSquares[i] = avgSquares[i] + (point[i]*point[i]-avgSquares[i])/newWeight;
    }
  }

  editTime = timestamp;
  weight = newWeight;
}

std::vector<double> microCluster::getCenter(){
  return center;
}

int microCluster::getCreationTime(){
  return creationTime;
}

double microCluster::getWeight(int timestamp){ //Weight decays depending on access time
  return weight * pow(2, -lambda*(timestamp-editTime));
}

double microCluster::getRadius(){
  double r2 = 0;
  for(int i=0; i<center.size(); i++){
    r2 += avgSquares[i] - center[i]*center[i]; //(CF2/w)-(CF1/w)^2
  }

  return sqrt(r2);
}

double microCluster::mergeRadius(std::vector<double> point, int timestamp){ //Radius if new point would be merged in
  double newWeight = weight*pow(2, -lambda*(timestamp - editTime)) + 1;
  double r2 = 0;

  for(int i=0; i<center.size(); i++){
    double newCenter = center[i] + (point[i]-center[i])/newWeight;
    double newAvgSquare = avgSquares[i] + (point[i]*point[i]-avgSquares[i])/newWeight;
    r2 += newAvgSquare - newCenter*newCenter; //Use updated center and average sum of squares
  }

  return sqrt(r2);
}

// basePipe constructor
template <typename nodeType>
denStream<nodeType>::denStream(){
	this->procName = "DenStream";
    return;
}
//taking in preprocessor type


// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void denStream<nodeType>::runPreprocessor(pipePacket<nodeType>& inData){
  /////////constants//////////
  initPoints = 1000; // points to generate p clusters (large data sets, use 1000) 
  minPoints = 20; //For DBSCAN 
  mu = 10;   //weight of data points in cluster threshold
  beta = 0.2; // outlier threshold
  streamSpeed = 20; //number of points in one unit time

  Tp = ceil((1/lambda)*log((beta*mu)/((beta*mu)-1)));
  timestamp = 0; //Current time
  int numPerTime = 0; //Number of points processed so far in this unit time
  double xi_den = pow(2, -Tp*lambda)-1; //Denominator for outlier lower weight limit

  ///////initialize p micro clusters... DBSCAN first N points (N has to be less than size of input data to simulate stream)///////
  //this returns cluster labels corresponding to current points
  std::vector<int> clusterLabels = dbscan::cluster(inData.workData, minPoints, epsilon, initPoints);
  int pClusters = *std::max_element(clusterLabels.begin(), clusterLabels.end()); //Number of clusters found

  timestamp = initPoints/streamSpeed;
  numPerTime = initPoints % streamSpeed;

  if(pClusters == -1) pClusters = 1;
  pMicroClusters = std::vector<microCluster>(pClusters, microCluster(timestamp, inData.workData[0].size(), lambda));

  //Take labels of each point and insert into appropriate microCluster
  for(int i = 0; i<initPoints; i++){
    if(clusterLabels[i] != -1){ //Not noise
      pMicroClusters[clusterLabels[i]-1].insertPoint(inData.workData[i], timestamp);
    }
  }

  for(int i = initPoints; i<inData.workData.size(); i++){
    numPerTime++;
    if(numPerTime == streamSpeed){ //1 unit time has elapsed
      numPerTime = 0;
      timestamp++;
    }

    merging(inData.workData, i, timestamp);

    if(timestamp % Tp == 0 && numPerTime == 0){ //Check microclusters after Tp units of time
      auto itP = pMicroClusters.begin();
      while(itP != pMicroClusters.end()){
        if(itP->getWeight(timestamp) < beta*mu){ //Weight too low - no longer a pMicroCluster
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

        if(itO->getWeight(timestamp) < xi){ //Weight too low - no longer an oMicroCluster
          itO = oMicroClusters.erase(itO);
        } else{
          ++itO;
        }
      }
    }
  }

  std::vector<std::vector<double>> centers(pMicroClusters.size()); //Return p-micro-cluster centers
  for(int i=0; i<pMicroClusters.size(); i++){
    centers[i] = pMicroClusters[i].getCenter();
  }

  inData.workData = centers;
	return;
}

//Index of nearest p-micro-cluster
template <typename nodeType>
int denStream<nodeType>::nearestPCluster(std::vector<double> point){
  utils ut;

  int min = 0;
  double minDist = ut.vectors_distance(pMicroClusters[0].getCenter(), point); 

  for(int i=1; i<pMicroClusters.size(); i++){
    double dist = ut.vectors_distance(pMicroClusters[i].getCenter(), point); //Find microcluster with least distance from center to point
    if(dist < minDist){
      minDist = dist;
      min = i;
    }
  }

  return min;
}

//Index of nearest o-micro-cluster
template <typename nodeType>
int denStream<nodeType>::nearestOCluster(std::vector<double> point){
  utils ut;

  int min = 0;
  double minDist = ut.vectors_distance(oMicroClusters[0].getCenter(), point); 

  for(int i=1; i<oMicroClusters.size(); i++){
    double dist = ut.vectors_distance(oMicroClusters[i].getCenter(), point); //Find microcluster with least distance from center to point
    if(dist < minDist){
      minDist = dist;
      min = i;
    }
  }

  return min;
}

//function merging (takes in point p, p micro clusters, o micro clusters )
template <typename nodeType>
void denStream<nodeType>::merging(std::vector<std::vector<double>> &data, int p, int timestamp){
  if(pMicroClusters.size() > 0){
    int nearestP = nearestPCluster(data[p]);
    if(pMicroClusters[nearestP].mergeRadius(data[p], timestamp) <= epsilon){ //Successfully add point to nearest p-micro-cluster
      pMicroClusters[nearestP].insertPoint(data[p], timestamp);
      return;
    }
  }

  if(oMicroClusters.size() > 0){
    int nearestO = nearestOCluster(data[p]);
    if(oMicroClusters[nearestO].mergeRadius(data[p], timestamp) <= epsilon){ //Successfully add point to nearest o-micro-cluster
      oMicroClusters[nearestO].insertPoint(data[p], timestamp);
      if(oMicroClusters[nearestO].getWeight(timestamp) > beta*mu){ //Move from o-micro-cluster to p-micro-cluster
        pMicroClusters.push_back(oMicroClusters[nearestO]);
        oMicroClusters.erase(oMicroClusters.begin()+nearestO);
      }
      return;
    }
  }

  microCluster x(timestamp, data[0].size(), lambda); //Create new o-micro-cluster
  x.insertPoint(data[p], timestamp);
  oMicroClusters.push_back(x);
}
  

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool denStream<nodeType>::configPreprocessor(std::map<std::string, std::string> &configMap){
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
  
  pipe = configMap.find("epsilon"); //Maximum radius of micro clusters
  if(pipe !=configMap.end())
    this->epsilon = std::stod(configMap["epsilon"].c_str());
  else return false;

  pipe = configMap.find("lambda"); //Decay factor
  if(pipe !=configMap.end())
    this->lambda = std::stod(configMap["lambda"].c_str());
  else return false;

  this->configured = true;
  this->ut.writeDebug("DenStream","Configured with parameters { epsilon: " + std::to_string(this->epsilon) + ", debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

template class denStream<simplexNode>;
template class denStream<alphaNode>;
template class denStream<witnessNode>;
