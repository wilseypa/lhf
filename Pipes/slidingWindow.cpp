/*
 * streamVR hpp + cpp extend the basePipe class for
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include <functional>
#include "slidingWindow.hpp"
#include "utils.hpp"
#include "readInput.hpp"


// basePipe constructor
slidingWindow::slidingWindow(){
	pipeType = "SlidingWindow";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket slidingWindow::runPipe(pipePacket inData){
	utils ut;
	readInput rp;

	int windowMaxSize = 200;
	int windowMinSize = 150;

	double minNNdist = 0.5;
	double maxNNdist = 5.0;

	// For this pipe, we construct a sub-pipeline:
	//		1. Read data vector by vector, push into slidingWindow evaluation
	//		2. IF the point is to be inserted, push into the simplexTree (LRU ordered)
	//		3. IF 100 points have passed, generate persistence intervals
	//
	// Loop this subpipeline until there's no more data

	std::vector<int> windowKeys;
	std::vector<std::vector<double>> windowValues;
	std::vector<int> dynamicKeyContainer;

	std::vector<double> currentVector;
	rp.streamInit(inputFile);
	int pointCounter = 1;
	int key = 0;
	int indexCounter = 0;

	while(rp.streamRead(currentVector)){

		//Evaluate insertion into sliding window
		if(windowKeys.size() < windowMinSize){
			windowKeys.push_back(key);
			windowValues.push_back(currentVector);
			dynamicKeyContainer.push_back(key);
			key++;

			//If we've reached window size, generate the initial complex
			if(key == windowMinSize){
				std::cout << "Initializing complex" << std::endl;

				inData.originalData = windowValues;
				runComplexInitializer(inData);

				std::cout << "Returning from complex initializer" << std::endl;
			}

		} else {
			//Do NN search

			//One with LHF aNN (Aaron)
			//  auto [index, sqrdDistance] = ut.aNN(dataPoint, windowValues);

			//One with ANN aNN (Streaming branch)
			//	auto [index, sqrdDistance] = nnSearch(dataPoint, windowValues);

			//One with LHF NN (Utils)
			auto reps = ut.nearestNeighbors(currentVector, windowValues);

			//ut.print1DVector(reps);

			double minDist = *std::min_element(reps.begin(),reps.end());
			double maxDist = *std::max_element(reps.begin(),reps.end());

			//std::cout << "\tMin: " << minDist << "\tMax: " << maxDist << "\tBounds ( " << minNNdist << " , " << maxNNdist << ")" << std::endl;

			double nearRep = minDist;

			if((nearRep < minNNdist) || (nearRep > maxNNdist)){
				indexCounter++;

				//Check if we need to remove points
				if(windowKeys.size() == windowMaxSize) {
					std::cout << "\tDeleting..." << std::endl;

					int keyToDelete = dynamicKeyContainer[0];
					dynamicKeyContainer.erase(dynamicKeyContainer.begin());

					auto it = std::find(windowKeys.begin(), windowKeys.end(), keyToDelete);
					int indexToDelete = std::distance(windowKeys.begin(), it);

					windowKeys.erase(windowKeys.begin() + indexToDelete);
					windowValues.erase(windowValues.begin() + indexToDelete);
					inData.complex->deleteIterative(keyToDelete);

				}

				//Insert the point
				inData.complex->insertIterative(reps);
				windowValues.push_back(currentVector);
				windowKeys.push_back(key);
				dynamicKeyContainer.push_back(key);
				key++;

				//Still need to update min/max NN distances
//				if(minDist < minNNdist)
//					minNNdist = minDist;
//				if(maxDist > maxNNdist)
//					maxNNdist = maxDist;

			} else {
				int nnIndex = 0;//index[0];
				int nnKey = windowKeys[(unsigned)nnIndex];

				auto itRep = std::find(dynamicKeyContainer.begin(), dynamicKeyContainer.end(), nnKey);

				std::rotate(itRep, itRep + 1, dynamicKeyContainer.end());

			}
		}

		//Check if we've gone through 100 points
		if(pointCounter % 100 == 0 && pointCounter > 100){
			// Build and trigger remaining pipeline. It should only require the computation of persistence
			// intervals from the complex being maintained.
			std::cout << "pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;

			runSubPipeline(inData);

		}
		pointCounter++;

		currentVector.clear();

	}
	//Probably want to trigger the remaining pipeline one last time...
	if((pointCounter - 1) % 100 != 0){
		std::cout << "pointCounter: " << pointCounter << "\tSimplex Count: " << inData.complex->simplexCount() << "\tVertex Count: " << inData.complex->vertexCount() << std::endl;

		runSubPipeline(inData);
	}


	ut.writeLog("slidingWindow", "\tSuccessfully evaluated " + std::to_string(pointCounter) + " points");

	return inData;
}


void slidingWindow::runSubPipeline(pipePacket wrData){
    if(wrData.originalData.size() == 0)
		return;

	pipePacket inData = wrData;

	std::string pipeFuncts = "rips.persistence";
	auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;

	//For each '.' separated pipeline function (count of '.' + 1 -> lim)
	for(unsigned i = 0; i < lim; i++){
		auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
		pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

		//Build the pipe component, configure and run
		auto *bp = new basePipe();
		auto *cp = bp->newPipe(curFunct, "simplexTree");

		//Check if the pipe was created and configure
		if(cp != 0 && cp->configPipe(subConfigMap)){
			//Run the pipe function (wrapper)
			inData = cp->runPipeWrapper(inData);
		} else {
			std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
		}
	}

	return;
}

void slidingWindow::runComplexInitializer(pipePacket &inData){
	if(inData.originalData.size() == 0)
		return;

	std::string pipeFuncts = "distMatrix.neighGraph";
	auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;
	//For each '.' separated pipeline function (count of '.' + 1 -> lim)
	for(unsigned i = 0; i < lim; i++){
		auto curFunct = pipeFuncts.substr(0,pipeFuncts.find('.'));
		pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

		//Build the pipe component, configure and run
		auto *bp = new basePipe();
		auto *cp = bp->newPipe(curFunct, "simplexTree");

		//Check if the pipe was created and configure
		if(cp != 0 && cp->configPipe(subConfigMap)){
			//Run the pipe function (wrapper)
			inData = cp->runPipeWrapper(inData);
		} else {
			std::cout << "LHF : Failed to configure pipeline: " << curFunct << std::endl;
		}
	}
	return;
}



// configPipe -> configure the function settings of this pipeline segment
bool slidingWindow::configPipe(std::map<std::string, std::string> configMap){
	std::string strDebug;
	subConfigMap = configMap;

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = configMap["outputFile"].c_str();

	ut = utils(strDebug, outputFile);

	pipe = configMap.find("inputFile");
	if(pipe != configMap.end())
		inputFile = configMap["inputFile"].c_str();

	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		epsilon = std::atof(configMap["epsilon"].c_str());
	else return false;

	pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	else return false;

	ut.writeDebug("slidingWindow","Configured with parameters { input: " + configMap["inputFile"] + ", dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");

	return true;
}


// outputData -> used for tracking each stage of the pipeline's data output without runtime
void slidingWindow::outputData(pipePacket inData){
	std::ofstream file ("output/" + pipeType + "_output.csv");

	if(inData.complex->simplexType == "simplexArrayList"){
		for(auto a : inData.complex->weightedGraph[1]){
			for(auto d : a.first){
				file << d << ",";
			}
			file << "\n";
		}
	}else{
		auto edges = inData.complex->getAllEdges(5);
		for (auto edge : edges){
			for (auto a : edge){
				for(auto d : a.first){
					file << d << ",";
				}
				file << a.second << "\n";
			}
		}
	}
	file << std::endl;

	file.close();
	return;
}
