#include "directPolytopal.hpp"
#include "utils.hpp"
#include <bits/stdc++.h>
#include <iostream>
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/PointCoordinates.h"
// basePipe constructor
template <typename nodeType>
directPolytopal<nodeType>::directPolytopal()
{
  this->pipeType = "directPolytopal";
  return;
}
// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void directPolytopal<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
  double mini=DBL_MAX;
  double maxi=DBL_MIN;

  std::vector<std::vector <double>> data = inData.workData;
  for (auto x : data) {
      for (auto y : x) {
          std::cout << y << ",";
      }
      std::cout << std::endl;
  }
  std::cout<<"Yogesh will write evolutionary algorithm here"<<std::endl;

  std::vector<std::pair<double, double>> range_min_max;


  for(int i=0;i<data.size();i++)
  {
    for(int j=0;j<data[i].size();j++)
    {
      mini=std::min(mini,data[i][j]);
      maxi=std::max(maxi,data[i][j]);
    }
    range_min_max.push_back({mini,maxi});
  }


  std::vector<std::vector<double>> random_points;
    int m;
    std::cout << "Enter the number of random points to insert: ";
    std::cin >> m;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist;

    for (int i = 0; i < m ; i++)
    {
      std::vector<double> rpoint;
      for(int j=0;j<data[0].size();j++)
      {
        dist = std::uniform_real_distribution<>(range_min_max[j].first, range_min_max[j].second);
        rpoint.push_back(dist(gen));
      }
      
      random_points.push_back(rpoint);
    }



    std::cout << "Generated points:\n";
    for (auto p : data)
    {
      // std::cout << "(" << p.first << ", " << p.second << ")\n";
      for(int i=0;i<p.size();i++)
      {
        std::cout << p[i] << ",";
      }
      std::cout << std::endl;
    }
    std::cout << "Random points to insert:\n";
    for (auto p : random_points)
    {
      for(int i=0;i<p.size();i++)
      {
        std::cout << p[i] << ",";
      }
      std::cout << std::endl;
    }

    // ------------------ FIND NEAREST POINT ------------------
    double radius = maxi;
    for (auto p : data)
    {
      double d = utils::vectors_distance(p, random_points[0]); 
      if (d < radius)
      {
        radius = d;
      }
    }

    // ------------------ INVERSION ------------------
    std::vector<std::vector<double>> inverted_points;

    for (auto p : data)
    {
      double dist_project_point = utils::vectors_distance(p, random_points[0]);

      double factor = (radius * radius) / (dist_project_point * dist_project_point);
      std::vector<double> proj_point;
      for(int i=0;i<p.size();i++)      {
        proj_point.push_back(random_points[0][i] + factor * (p[i] - random_points[0][i]));
      } 
      
      inverted_points.push_back(proj_point);
    }

    std::cout<<"Radius of inversion: "<<radius<<std::endl; 
    std::cout<<"random point: ";
    for(int i=0;i<random_points[0].size();i++)    {
      std::cout << random_points[0][i] << ",";
    }
    std::cout << std::endl;

    std::cout<<"Inverted points:\n";

    for(auto p : inverted_points)
    {
      for(int i=0;i<p.size();i++)
      {
        std::cout << p[i] << ",";
      }
      std::cout << std::endl;
    }

    // ------------------ CONVEX HULL OF INVERTED POINTS ------------------
    std::cout << "Computing Convex Hull..." << std::endl;
    orgQhull::Qhull qh;
    std::vector<double> sdata;
    for (auto& pt : inverted_points)
        for (auto& coord : pt)
            sdata.push_back(coord);

    orgQhull::PointCoordinates pts(qh, data[0].size(), "Inverted Points");
    pts.append(sdata);
    std::cout<<"  pts.comment(): "<<pts.comment().c_str()<<std::endl;
    std::cout<<"  pts.dimension(): "<<pts.dimension()<<std::endl;
    std::cout<<"  pts.count(): "<<pts.count()<<std::endl;

    qh.runQhull(pts.comment().c_str(), pts.dimension(), pts.count(),&*pts.coordinates(), "Qt");

    // Extract hull boundary vertices and facets
    std::set<unsigned> hull_vertex_ids;
    std::vector<std::vector<unsigned>> hull_facets;

    for (orgQhull::QhullFacet f : qh.facetList()) {
        std::vector<unsigned> facet_vertices;
        for (orgQhull::QhullVertex v : f.vertices()) {
            unsigned id = v.point().id();
            hull_vertex_ids.insert(id);
            facet_vertices.push_back(id);
        }
        hull_facets.push_back(facet_vertices);
    }
    // Print hull results
    std::cout << "Convex Hull: " << hull_vertex_ids.size()<< " vertices, " << hull_facets.size() << " facets" << std::endl;
    std::cout << "Hull vertex indices: ";
    for (auto id : hull_vertex_ids)
        std::cout << id << " ";
    std::cout << std::endl;

    std::cout << "Hull boundary points (original data):" << std::endl;
    for (auto id : hull_vertex_ids) {
        for (int i = 0; i < data[id].size(); i++)
            std::cout << data[id][i] << ",";
        std::cout << std::endl;
    }

  /*Tentative algorithm
 A) Generate initial population using following procedure to generate all N of them:

  1. Generate (r) random points within in a point cloud
  2. For each r_i in r:
  3.    Find k-nearest neighbors in the vicinity of a random point r_i
  4.    Find the nearest point(min distance r_min) in the K-NN to r_i and assume a hypersphere of radius r_min
  5.    inverse Project each point(inverse_i) from knn into this hypersphere using parametric equation x_p = r_i(1-r_min^2/d_i^2)+ inverse_i(r_min^2/d_i^2)
  6.    Find the approximate outer boundary, and keep that as a hull boundary on original point cloud
  7. Keep the maximal non overlapping convex polytopes
 B) Find the fitness of each member of the population based on number of points covered and maximal size od the polytopes

 C) Generate new population using crossover and mutation that has better fitness

 D) Keep the fitest offsprings for next iteration repeate C and D until we find the best result.

  */
  return;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool directPolytopal<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
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

  this->ut = utils(strDebug, this->outputFile);

  this->configured = true;
  this->ut.writeDebug("directPolytopal", "Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

  return true;
}
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void directPolytopal<nodeType>::outputData(pipePacket<nodeType>& inData)
{
    // Output related to betaSubSkeletonComplex
    return;
}
template class directPolytopal<simplexNode>;
template class directPolytopal<alphaNode>;
template class directPolytopal<witnessNode>;
