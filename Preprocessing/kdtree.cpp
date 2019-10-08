*
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
#include "kdtree.hpp"
#include "utils.hpp"

// basePipe constructor
kdtree::kdtree()
{
    procName = "kdtree";
    return;
}

// runPipe -> Run the configured functions of this pipeline segment
pipePacket kdtree::runPreprocessor(pipePacket inData)
{

}

// configPipe -> configure the function settings of this pipeline segment
bool kdtree::configPreprocessor(std::map<std::string, std::string> configMap)
{
    std::string strDebug;

    auto pipe = configMap.find("debug");
    if (pipe != configMap.end())
    {
        debug = std::atoi(configMap["debug"].c_str());
        strDebug = configMap["debug"];
    }
    pipe = configMap.find("outputFile");
    if (pipe != configMap.end())
        outputFile = configMap["outputFile"].c_str();


  //  ut.writeDebug("StreamKMeans", "Configured with parameters { clusters: " + configMap["clusters"] + ", iterations: " + configMap["iterations"] + ", debug: " + strDebug + ", outputFile: " + outputFile + " }");

}
