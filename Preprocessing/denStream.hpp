#pragma once

// Header file for denStream class - see denStream.cpp for descriptions
#include <map>
#include <vector>
#include "preprocessor.hpp"

class microCluster{
	private:
		std::vector<double> center;
		std::vector<double> avgSquares;
		double weight;
		double lambda;
		int creationTime;
		int editTime;

	public:
		microCluster(int t, int dim, double lambda);
		void insertPoint(std::vector<double> point, int timestamp);
		std::vector<double> getCenter();
		int getCreationTime();
		double getRadius(int timestamp);
		double getWeight(int timestamp);
		double mergeRadius(std::vector<double> point, int timestamp);
};

class denStream : public preprocessor {
 	private:
		int initPoints;
		int minPoints;
		double lambda;
		double epsilon;
		int mu;
		double beta;    
		int streamSpeed;
		int Tp;
		int timestamp;

		std::vector<microCluster> pMicroClusters;
  		std::vector<microCluster> oMicroClusters;

		void merging(std::vector<std::vector<double>> &data, int p, int timestamp);
		int nearestPCluster(std::vector<double> point);
		int nearestOCluster(std::vector<double> point);
  	
	public:
		denStream();
		pipePacket runPreprocessor(pipePacket inData);
		bool configPreprocessor(std::map<std::string, std::string> configMap);
}; 
