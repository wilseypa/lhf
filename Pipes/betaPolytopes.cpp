/*
 * betaSubSkeletonComplexPipe hpp + cpp extend the basePipe class for calculating the
 * beta Skeleton Based Complex generation for data input
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <set>
#include <algorithm>
#include "betaPolytopes.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include "readInput.hpp"

// basePipe constructor
template <typename nodeType>
betaPolytopes<nodeType>::betaPolytopes()
{
	this->pipeType = "betaPolytopes";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void betaPolytopes<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool betaPolytopes<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
{
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if (pipe != configMap.end())
	{
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if (pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();

	pipe = configMap.find("beta");
	if (pipe != configMap.end())
		this->beta = std::atof(configMap["beta"].c_str());

	pipe = configMap.find("betaMode");
	if (pipe != configMap.end())
		this->betaMode = configMap["betaMode"].c_str();

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->epsilon = std::atof(configMap["epsilon"].c_str());

	this->ut = utils(strDebug, this->outputFile);
	pipe = configMap.find("dimensions");
	if (pipe != configMap.end())
	{
		this->dim = std::atoi(configMap["dimensions"].c_str());
	}
	pipe = configMap.find("betaMesh");
	if (pipe != configMap.end())
		this->betaMesh = configMap["betaMesh"].c_str();

	pipe = configMap.find("epsilon");
	if (pipe != configMap.end())
		this->enclosingRadius = std::atof(configMap["epsilon"].c_str());
	else
		return false;

	this->configured = true;
	this->ut.writeDebug("betaSubSkeletonComplex Pipe ", "Configured with parameters { eps: " + configMap["epsilon"] + configMap["beta"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}

// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void betaPolytopes<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	// Output related to betaPolytopes
	return;
}


template class betaPolytopes<simplexNode>;
template class betaPolytopes<alphaNode>;
template class betaPolytopes<witnessNode>;
