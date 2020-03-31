#pragma once

// Header file for kdTree class - see kdTree.cpp for descriptions
#include <map>
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

        double coord(const size_t &);
        
        explicit operator bool();
        explicit operator point();
        explicit operator size_t();
        explicit operator pointIndex();
};

using kdNodePtr = std::shared_ptr< kdNode >;

kdNodePtr newkdNodePtr();

using pointIndexArr = typename std::vector<pointIndex>;
using pointVec = std::vector<point>;

class kdTree{
    public: 
        kdNodePtr root;

        kdNodePtr makeTree(const pointIndexArr::iterator &begin, const pointIndexArr::iterator &end, 
                       const size_t &length, const size_t &level);
    public:
        kdTree();
        //explicit kdTree(pipePacket inData); //Prevent implicit conversion
        kdTree(pointVec inData, int size);

    private:
        kdNodePtr findNearest(
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
        pointIndexArr neighborhood(
            const kdNodePtr &branch,
            const point &pt,
            const double &rad, //epsilon
            const size_t &level
        );

    public:
        pointIndexArr neighborhood(const point &pt, const double &rad);

        pointVec neighborhoodPoints(const point &pt, const double &rad);

        indexArr neighborhoodIndices(const point &pt, const double &rad);
};
