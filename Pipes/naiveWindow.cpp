/*
 * naiveWindow hpp + cpp extend the basePipe class for
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include <cmath>
#include <numeric>
#include <functional>
#include "naiveWindow.hpp"
#include "readInput.hpp"

// basePipe constructor
template <typename nodeType>
naiveWindow<nodeType>::naiveWindow()
{
	this->pipeType = "NaiveWindow";
	return;
}

template <typename nodeType>
bool naiveWindow<nodeType>::sampleStreamEvaluator(std::vector<double> &vector, std::vector<std::vector<double>> &window)
{
	// Return true to accept all points into the complex
	return true;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void naiveWindow<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	readInput rp;

	int windowMaxSize = 100;

	// For this pipe, we construct a sub-pipeline:
	//		1. Read data vector by vector, push into stream evaluation in complex
	//		2. IF the point is to be inserted, push into the simplexTree (LRU ordered)
	//		3. IF 100 points have passed, generate persistence intervals
	//
	// Loop this subpipeline until there's no more data

	int key = 0;
	std::vector<int> windowKeys;
	std::vector<std::vector<double>> windowValues;
	std::vector<double> currentVector;
	int indexToBeDeleted = 0;

	// Initialize the read stream using the input file
	if (rp.streamInit(inputFile))
	{

		int pointCounter = 1;
		// Get the next point from the stream
		while (rp.streamRead(currentVector))
		{
			// Evaluate insertion into sliding window

			// Window hasn't filled
			if (windowValues.size() < windowMaxSize)
			{
				windowValues.push_back(currentVector);
				windowKeys.push_back(key);
				key++;
				// If we've reached window size, generate the initial complex
				if (windowValues.size() == windowMaxSize)
				{
					inData.workData = windowValues;

					distMatrix.resize(inData.workData.size(), std::vector<double>(inData.workData.size(), 0));

					// Iterate through each vector
					for (unsigned i = 0; i < inData.workData.size(); i++)
					{
						if (!inData.workData[i].empty())
						{
							// Grab a second vector to compare to
							std::vector<double> temp;
							for (unsigned j = i + 1; j < inData.workData.size(); j++)
							{
								// Calculate vector distance
								auto dist = this->ut.vectors_distance(inData.workData[i], inData.workData[j]);

								if (dist < epsilon)
									distMatrix[i][j] = dist;
							}
						}
					}

					// Set the stream evaluator
					// inData.complex->setStreamEvaluator(&this->sampleStreamEvaluator);
				}

				// Window is full, evaluate and add to window
			}
			else
			{
				if (inData.complex->insertIterative(currentVector, windowValues))
				{

					windowValues.erase(windowValues.begin());

					// Insert the point
					windowValues.push_back(currentVector);
					windowKeys.push_back(key);
					key++;
				}
			}

			// Check if we've gone through 100 points
			if (pointCounter % 10 == 0 && pointCounter >= windowMaxSize)
			{
				// Build and trigger remaining pipeline. It should only require the computation of persistence
				// intervals from the complex being maintained.
				inData.workData = windowValues;
				runSubPipeline(inData);
				inData = pipePacket<nodeType>("simplexArrayList", epsilon, dim);
			}
			pointCounter++;

			currentVector.clear();
		}
		// Probably want to trigger the remaining pipeline one last time...
		if ((pointCounter - 1) % 10 != 0)
		{
			runSubPipeline(inData);
		}
		this->ut.writeLog("naiveWindow", "\tSuccessfully evaluated " + std::to_string(pointCounter) + " points");

		writeComplexStats(inData);
	}

	return;
}

template <typename nodeType>
void naiveWindow<nodeType>::writeComplexStats(pipePacket<nodeType> &inData)
{
	if (inData.complex->stats.size() > 30)
	{
		std::ofstream file("output/complexStats.csv");

		file << inData.complex->stats << std::endl;

		file.close();
	}
}

template <typename nodeType>
void naiveWindow<nodeType>::runSubPipeline(pipePacket<nodeType> wrData)
{
	if (wrData.workData.size() == 0)
		return;

	pipePacket<nodeType> inData = wrData;
	outputData(inData);

	inData.complex->setDistanceMatrix(&distMatrix);
	std::string pipeFuncts = "neighGraph.incrementalPersistence";
	subConfigMap["complexType"] = "simplexArrayList";
	auto lim = count(pipeFuncts.begin(), pipeFuncts.end(), '.') + 1;
	subConfigMap["fn"] = "_" + std::to_string(repCounter);

	this->repCounter++;

	// For each '.' separated pipeline function (count of '.' + 1 -> lim)
	for (unsigned i = 0; i < lim; i++)
	{
		auto curFunct = pipeFuncts.substr(0, pipeFuncts.find('.'));
		pipeFuncts = pipeFuncts.substr(pipeFuncts.find('.') + 1);

		// Build the pipe component, configure and run
		auto cp = basePipe<nodeType>::newPipe(curFunct, "simplexTree");

		// Check if the pipe was created and configure
		if (cp != 0 && cp->configPipe(subConfigMap))
		{
			// Run the pipe function (wrapper)
			cp->runPipeWrapper(inData);
		}
		else
		{
			std::cout << "LHF subPipe: Failed to configure pipeline: " << curFunct << std::endl;
		}
	}

	distMatrix.clear();

	return;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool naiveWindow<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
{
	std::string strDebug;
	subConfigMap = configMap;

	auto pipe = configMap.find("debug");
	if (pipe != configMap.end())
	{
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if (pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();

	this->ut = utils(strDebug, this->outputFile);

	pipe = configMap.find("inputFile");
	if (pipe != configMap.end())
		this->inputFile = configMap["inputFile"].c_str();

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->epsilon = std::atof(configMap["epsilon"].c_str());
	else
		return false;

	pipe = configMap.find("dimensions");
	if (pipe != configMap.end())
	{
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	else
		return false;

	this->configured = true;
	this->ut.writeDebug("naiveWindow", "Configured with parameters { input: " + configMap["inputFile"] + ", dim: " + configMap["dimensions"] + ", eps: " + configMap["epsilon"] + ", debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void naiveWindow<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	std::ofstream file("output/" + this->pipeType + "_" + std::to_string(this->repCounter) + "_output.csv");

	for (auto a : inData.workData)
	{
		for (auto d : a)
		{
			file << d << ",";
		}
		file << "\n";
	}
	file << std::endl;

	file.close();
	return;
}

template class naiveWindow<simplexNode>;
template class naiveWindow<alphaNode>;
template class naiveWindow<witnessNode>;
