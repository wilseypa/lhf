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
#include "kdTree.hpp"
#include "utils.hpp"

/*kdTree::kdTree(){
    procName = "kdTree";
    return;
} */
// create a binary tree of nodes of "k" dimensions from points of the input csv
//every leaf node is a k dimensional point, non leaf nodes splits hyperplane of space
// each node represents axis aligned hyper rectangle 
// each node specifies an axis and splits the set of points based on whether
//their coordinate along that axis is > || < a certain value
//axis and splitting point chosen by sliding midpoint rule -> split along longest side of node
// ie divide perpendicular to longest spread 
//docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.html

//kdtree is used to find nearest neighbors to current point based on epsilon




//input: List of points (inData), depth ( O log(n))
std::vector<double> a{ 1.1, 2.2, 1.7, 5.6, 1.123, 9.1};
std::vector<double> b{3.34, 5.23, .532, 1.12, 3.45, 1.23};
std::vector<std::vector<double>> pointList { a , b};
int depth = log2(pointList.size());
int k = pointList[0].size();

//select axis based on depth so that it cycles through all valid vals
int axis = depth%k; 


// sort pointList and choose median as pivot element
//std::sort(pointList[axis].begin(), pointList[axis].end());


//create node and construct subtree
/* node location = median
leftChild = kdtree( points in pointList before median, depth + 1)
rightChild = kdtree(points in pointList after median, depth +1)
return node en.wikipedia.org/wiki/K-d_tree#Construction
*/

//node has 5 fields: splitting axis, split val, L/R subtree, point (if leaf)






// runPipe -> Run the configured functions of this pipeline segment
pipePacket kdTree::runPreprocessor(pipePacket inData)
{

}




// configPipe -> configure the function settings of this pipeline segment
bool kdTree::configPreprocessor(std::map<std::string, std::string> configMap)
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
