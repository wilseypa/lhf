#pragma once

// Header file for kdTree class - see kdTree.cpp for descriptions
#include <map>
#include "preprocessor.hpp"
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

//github.com/crvs/KDTree/blob/master/KDTree.hpp
using point = std::vector<double>;
using pointList = std::vector<std::vector<double>>;
using pointListItr = pointList::iterator; 
using indexArr = std::vector< size_t >;
using pointIndex = typename std::pair<std::vector<double>, size_t>;  //alias declarations


class kdNode {
    public: 
        using kdNodePtr = std::shared_ptr< kdNode >; //initializer smart pointer
        size_t index;
        point x;
        kdNodePtr left; 
        kdNodePtr right;


        kdNode();
        kdNode(const point &, const size_t &, const kdNodePtr &, const kdNodePtr &);
        kdNode(const pointIndex &, const kdNodePtr &, const kdNodePtr &);
        ~kdNode(); //destructor

};

using kdNodePtr = std::shared_ptr< kdNode >;

kdNodePtr newkdNodePtr();

using pointIndexArr = typename std::vector<pointIndex>;
using pointVec = std::vector<point>;

class kdTree{
    kdNodePtr root;
    kdNodePtr leaf;

    kdNodePtr makeTree(const pointIndexArr::iterator &begin, const pointIndexArr::iterator &end, 
                       const size_t &length, const size_t &level);
    public:
        kdTree() = default;
       explicit kdTree(pointVec pointArray, pipePacket inData); //vector of points, which are std::vector<doubles>. prevent implicit conversion

    private:
      kdNodePtr nearest(
          const kdNodePtr &branch,
          const point &pt,
          const size_t &level,
          const kdNodePtr &best,
          const double &bestDist
      );

    kdNodePtr nearest(const point &pt);

    public:
        point nearestPoint(const point &pt);
        size_t nearestIndex(const point &pt);
        pointIndex nearestPointIndex(const point &pt);
    
    private:
        pointIndexArr neighborhood_(
            const kdNodePtr &branch,
            const point &pt,
            const double &rad, //epsilon
            const size_t &level
        );

    public:
        pointIndexArr neighborhood( //
            const point &pt,      //
            const double &rad);

        pointVec neighborhood_points( //
            const point &pt,        //
            const double &rad);

        indexArr neighborhood_indices( //
            const point &pt,         //
            const double &rad);
};

/*class kdTree : public preprocessor
{
private:
  

public:
    void makePointList(pipePacket inData);
    kdNodePtr makeTree(const pointIndexArr::iterator &begin, //
                       const pointIndexArr::iterator &end,   //
                       const size_t &length,                 //
                       const size_t &level);

    pipePacket runPreprocessor(pipePacket inData);
    bool configPreprocessor(std::map<std::string, std::string> configMap);
};  */



