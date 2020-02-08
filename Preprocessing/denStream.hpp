#pragma once

// Header file for denStream class - see denStream.cpp for descriptions
#include <map>
#include <vector>
#include "preprocessor.hpp"

class microCluster{
	private:
		std::vector<double> CF1;
		std::vector<double> CF2;
		double weight;
		int creationTime;

	public:
		microCluster(int t);
		void insertPoint(std::vector<double> point);
		std::vector<double> getCenter();
		double getRadius();
		double getWeight();
		double mergeRadius(std::vector<double>);
		int getCreationTime();
};

class denStream : public preprocessor {
 	private:
		int initPoints = 20;
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

		void merging(std::vector<std::vector<double>> &data, int p);
		int nearestPCluster(std::vector<std::vector<double>> &data, int p);
		int nearestOCluster(std::vector<std::vector<double>> &data, int p);
  	
	public:
		denStream();
		pipePacket runPreprocessor(pipePacket inData);
		bool configPreprocessor(std::map<std::string, std::string> configMap);
}; 
