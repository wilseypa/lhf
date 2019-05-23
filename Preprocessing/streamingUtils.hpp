#pragma once

// Header file for streamingUtils class - see streamingkMeans.cpp for descriptions
#include <map>
#include "preprocessor.hpp"


class streamingUtils : public preprocessor {
  struct Point;

  struct Facility {
      Point * p;
      float * totals;

   //   Point * samples[Max_Samples];
      int numSamples;
  };

  void init(Facility * f, int dim);   // make facility for point with dimension dim

  void copyFacility(Facility * a, Facility * b);  //deep copy a= b

  Facility * findNearest(Point * p, Facility *fs[], int n, float & distSQ);
  // return point to nearest facility to p from given list
  // sqaure of dist between 2 points and nearest facility returned by distSQ

 void assign(Point * p, Facility * f);  //assign point p to facility f

 void assign(Facility *a, Facility * b);  //assigns facility b to facility a

 void setToCoM(Facility * f);  
 // set facility designation point to be CoM of all points in that facility

void ballKmeans(Facility *f, double dist);
// ball kmeans

int findNN(Facility *fs[], int n, int i);
//find NN of fs[i]

};