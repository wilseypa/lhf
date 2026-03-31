#include "directPolytopal.hpp"
#include "utils.hpp"

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
  std::vector<std::vector <double>> data = inData.workData;
  for (auto x : data) {
      for (auto y : x) {
          std::cout << y << ",";
      }
      std::cout << std::endl;
  }
  std::cout<<"Yogesh will write evolutionary algorithm here"<<std::endl;

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
