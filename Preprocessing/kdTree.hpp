#pragma once

// Header file for kdTree class - see kdTree.cpp for descriptions
#include <map>
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

// github.com/crvs/KDTree/blob/master/KDTree.hpp
// alias declarations

class kdNode
{
public:
    size_t index;
    std::vector<double> x;
    std::shared_ptr<kdNode> left;
    std::shared_ptr<kdNode> right;

    kdNode();
    kdNode(const std::vector<double> &, const size_t &, const std::shared_ptr<kdNode> &, const std::shared_ptr<kdNode> &);
    kdNode(const std::pair<std::vector<double>, size_t> &, const std::shared_ptr<kdNode> &, const std::shared_ptr<kdNode> &);
    ~kdNode(); // destructor

    double coord(const size_t &);

    explicit operator bool();
    explicit operator std::vector<double>();
    explicit operator size_t();
    explicit operator std::pair<std::vector<double>, size_t>();
};

std::shared_ptr<kdNode> newkdNodePtr();

class kdTree
{
public:
    std::shared_ptr<kdNode> root;
    std::shared_ptr<kdNode> makeTree(const std::vector<std::pair<std::vector<double>, size_t>>::iterator &begin, const std::vector<std::pair<std::vector<double>, size_t>>::iterator &end, const size_t &length, const size_t &level);

    kdTree();
    // explicit kdTree(pipePacket inData); //Prevent implicit conversion
    kdTree(std::vector<std::vector<double>> inData, int size);
    std::vector<double> nearestPoint(const std::vector<double> &pt);
    size_t nearestIndex(const std::vector<double> &pt);
    std::pair<std::vector<double>, size_t> nearestPointIndex(const std::vector<double> &pt);
    std::vector<std::pair<std::vector<double>, size_t>> neighborhood(const std::vector<double> &pt, const double &rad);
    std::vector<std::vector<double>> neighborhoodPoints(const std::vector<double> &pt, const double &rad);
    std::vector<size_t> neighborhoodIndices(const std::vector<double> &pt, const double &rad);
    std::vector<std::vector<bool>> betaNeighbors(std::vector<std::vector<double>> &, double beta, std::string betaMode);

private:
    std::shared_ptr<kdNode> findNearest(const std::shared_ptr<kdNode> &branch, const std::vector<double> &pt, const size_t &level, const std::shared_ptr<kdNode> &best, const double &bestDist);
    std::shared_ptr<kdNode> nearest(const std::vector<double> &pt);
    std::vector<std::pair<std::vector<double>, size_t>> neighborhood(const std::shared_ptr<kdNode> &branch, const std::vector<double> &pt, const double &rad, const size_t &level);
};

