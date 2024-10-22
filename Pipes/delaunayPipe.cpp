#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "delaunayPipe.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include <fstream>
#include <omp.h>

// basePipe constructor
template <typename nodeType>
delaunayPipe<nodeType>::delaunayPipe()
{
  this->pipeType = "delaunayPipe";
  return;
}

std::vector<std::vector<unsigned>> qdelaunay_o(std::vector<std::vector<double>> &inputData)
{
  int dim = inputData[0].size();
  // typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dimension_tag<7>>> T;
  typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>> T;
  T dt(dim);
  std::vector<T::Point> points;
  unsigned index = 0;
  std::map<T::Point, unsigned> index_of_vertex;
  for (auto& i : inputData)
  {
    T::Point p(i.begin(), i.end());
    points.push_back(p);
    index_of_vertex[p] = index++;
  }
  dt.insert(points.begin(), points.end());
  points.clear();
  size_t number_of_finite_full_cell = dt.number_of_finite_full_cells();
  std::vector<std::vector<unsigned>> dsimplexes(number_of_finite_full_cell, std::vector<unsigned>(dim + 1));
#pragma omp parallel for
  for (int j = 0; j <= dim; ++j)
  {
    auto i = dt.finite_full_cells_begin();
    for (size_t k = 0; k < number_of_finite_full_cell; k++)
      dsimplexes[k][j] = index_of_vertex[(i++)->vertex(j)->point()];
  }
  return dsimplexes;
}
// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void delaunayPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
  std::vector<std::vector<unsigned>> dsimplexes = qdelaunay_o(inData.inputData);
  ((alphaComplex<alphaNode> *)inData.complex)->dsimplexmesh = dsimplexes;
  ((alphaComplex<nodeType> *)inData.complex)->buildAlphaComplex(dsimplexes, inData.inputData.size(), inData.inputData);
  return;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool delaunayPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
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
  this->ut.writeDebug("delaunayPipe", "Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

  return true;
}
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void delaunayPipe<nodeType>::outputData(pipePacket<nodeType> &inData)
{
  std::ofstream file;
  file.open("output/" + this->pipeType + "_output.csv");
  // code to print the data

  file.close();
  return;
}

template class delaunayPipe<simplexNode>;
template class delaunayPipe<alphaNode>;
template class delaunayPipe<witnessNode>;