#pragma once

// Header file for coresetUtils class - see streamingkMeans.cpp for descriptions
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <map>
#include "preprocessor.hpp"


class coresetUtils : public preprocessor {

   //////// point -> pipePacket point with additional info for coreset functionality
   struct point {  // metadata associated with pipe data in the coreset
       int dimension;

       double weight;  //clustering features
       double squareSum;
       float *coordinates;

       double curCost;
       int centerCost;

   };
   void initPoint(struct point *point, int dimension);  //initialize point
   void freePoint(struct point *point);  //delete point
   void copyPoint(struct point *firstPoint, struct point *secondPoint); //data of 1st pt to addr of 2nd
   void copyPointWithoutInit(struct point *firstPoint, struct point *secondPoint); //data of 1st pt to existing addr of 2nd

  //////////streaming coreset data structures
  struct bucket { //bucket holds a coreset of "m" points
    int curSize;
    struct point *points;
    struct point *spillover;  //describes points when a bucket reaches size m
  };

  struct bucketManager {  //metadata of buckets
    int numBuckets;
    int maxBucketSize;
    struct bucket *buckets;
  };

  //initialize tracker, n points, d dimensions
  void initTracker(struct bucketManager *manager, int n, int d, int maxSize);

  //insert point into bucket
  void insertPoint(struct point *p, struct bucketManager *manager);
  // if bucket is full extract coreset
  struct point  * getCoresetFromManager(struct bucketManager *manager, int d); };


/////////// coreset tree data structures
struct treeNode {
  int numPoints;

  struct point **points;
  struct point *center; //points to center of tree node
  struct treeNode *lc; //point to left child node
  struct treeNode *rc; //point to right child node
  struct treeNode *parent;  //point to parent node
  double cost;  //cost of treeNode for coreset extraction
};

//unions coresets A & B to make a new coreset of size k
void unionTreeCoreset(int k, int n1,int n2,int d, struct point *setA,struct point *setB,struct point *centers);