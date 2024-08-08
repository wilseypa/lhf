#pragma once

// Header file for unionFind class - see unionFind.cpp for descriptions
#include <vector>
/**
 * @brief
 *
 */
class unionFind
{
private:
    std::vector<int> rank, parent;

public:
    unionFind(int n);
    int find(int i);
    bool join(int x, int y);
};
