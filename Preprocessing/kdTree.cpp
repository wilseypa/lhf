/*
 * basePipe hpp + cpp protoype and define a base class for building
 * pipeline functions to execute
 *
 */
#include <stdio.h>
#include <fstream>
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
#include <queue>

#include "kdTree.hpp"
#include "utils.hpp"

kdTree::kdTree(){
 //   procName = "kdTree";
    return;
} 
// create a binary tree of nodes of "k" dimensions from points of the input csv
//every leaf node is a k dimensional point, non leaf nodes splits hyperplane of space
// each node represents axis aligned hyper rectangle 
// each node specifies an axis and splits the set of points based on whether
//their coordinate along that axis is > || < a certain value
//axis and splitting point chosen by sliding midpoint rule -> split along longest side of node
// ie divide perpendicular to longest spread 
//docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.html

//kdtree is used to find nearest neighbors to current point based on epsilon

//avdongre.wordpress.com/2011/06/14/kd-tree-in-c/
// github.com/crvs/KDTree/blob/master/KDTree.cpp#L41

using point =std::vector<double>;
using pointList = std::vector<std::vector<double>>;
using pointListItr =  pointList::iterator; //alias declarations

kdNode::kdNode() = default;

kdNode::kdNode(const point &pt, const size_t &idx_, const kdNodePtr &left_,
               const kdNodePtr &right_)
{
    x = pt;
    index = idx_;
    left = left_;
    right = right_;
}

kdNode::kdNode(const pointIndex &pi, const kdNodePtr &left_,
               const kdNodePtr &right_)
{
    x = pi.first;
    index = pi.second;
    left = left_;
    right = right_;
}

kdNode::~kdNode() = default; //destructor

kdNodePtr newKdNodePtr()
{
    kdNodePtr mynode = std::make_shared<kdNode>();
    return mynode;
}

kdNodePtr kdTree::makeTree(const pointIndexArr::iterator &begin, //
                           const pointIndexArr::iterator &end,   //
                           const size_t &length,                 //
                           const size_t &level){
    if(begin == end){
        return newKdNodePtr(); // sanity check for empty list
    }

    size_t dim = begin->first.size();

    if(length > 1){
        std::sort(begin, end, level);
    }

    auto middle = begin + (length/2);

    auto leftBegin = begin;
    auto leftEnd = middle;
    auto rightBegin = middle + 1;
    auto rightEnd = end;

    size_t leftLength = length/2;
    size_t rightLength = length - leftLength -1;

    kdNodePtr left; //recursively build left and right branches until leaf reached
    if (leftLength > 0 && dim > 0){
        left = makeTree(leftBegin, leftEnd, leftLength, (level+1) % dim);
    }
    else{
      //  left = leaf;
    }

    kdNodePtr right; 
    if(rightLength > 0 && dim > 0){
        right = makeTree(rightBegin, rightEnd, rightLength, (level+1)% dim);
    }
    else{
     //   right = leaf;
    }



    return std::make_shared< kdNode >(*middle, left, right);
    //returns pointer to root of the constructed kd tree




}

//driver code, push back inData.originalData to point vec so we can put in tree
kdTree::kdTree(pointVec pointArray, pipePacket inData){
    kdNodePtr leaf = std::make_shared< kdNode>();
    pointIndexArr arr;

    for(size_t i=0; i<inData.originalData.size(); i++){
        arr.push_back(pointIndex(inData.originalData[i],  i));
    }
    
    auto begin = arr.begin();
    auto end = arr.end();
    
    size_t length = arr.size();
    size_t level = 0; 

    root = kdTree::makeTree(begin, end, length, level);

}
























     /*pipePacket kdTree::runPreprocessor(pipePacket inData){

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

}  */
