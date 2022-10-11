#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include "delaunayPipe.hpp"
#include "utils.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <set>
#include <algorithm>
// basePipe constructor
template <typename nodeType>
delaunayPipe<nodeType>::delaunayPipe(){
	this->pipeType = "delaunayPipe";
	return;
}
// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void delaunayPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
  int dim = inData.inputData[0].size();
  typedef CGAL::Delaunay_triangulation<CGAL::Epick_d< CGAL::Dimension_tag<4> > >      T;
  T dt(dim);
  std::vector<T::Point> points;
  points.reserve(inData.inputData.size());
  for(auto i : inData.inputData){
    T::Point p(i.begin(), i.end());
    points.push_back(p);
  }
  T::Vertex_handle hint;
  int i = 0;
  for (std::vector<T::Point>::iterator it = points.begin(); it != points.end(); ++it) {
    if (T::Vertex_handle() != hint) {
      hint = dt.insert(*it, hint);
  }
    else {
      hint = dt.insert(*it);
    }
    printf("Processing: %d/%d\n", ++i, (int)points.size());
  }
  std::cout<<dt<<std::endl;
  /*
  std::cout<<"Number of vertices "<<dt.number_of_vertices()<<std::endl;
  for (auto i=dt.finite_vertices_begin();i!=dt.finite_vertices_end();i++)
    std::cout<<"Finite vertices "<<i->point()<<std::endl;
  std::cout<<"Number of finite full cells "<<dt.number_of_finite_full_cells()<<std::endl;
  for (auto i=dt.finite_full_cells_begin();i!=dt.finite_full_cells_end();i++){
    for (auto j=i->vertices_begin();j!=i->vertices_end();++j){
      std::cout<<j<<" ";}
    std::cout<<std::endl;}*/
  this->ut.writeDebug("delaunayPipe","\t Sucessfully Executed Pipe");
  return;
}
// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool delaunayPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap){
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end()){
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();

	this->ut = utils(strDebug, this->outputFile);

	this->configured = true;
	this->ut.writeDebug("delaunayPipe","Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void delaunayPipe<nodeType>::outputData(pipePacket<nodeType> &inData){
	std::ofstream file;
	file.open("output/" + this->pipeType + "_output.csv");
	//code to print the data
  
  
  
  file.close();
	return;
}

template class delaunayPipe<simplexNode>;
template class delaunayPipe<alphaNode>;
template class delaunayPipe<witnessNode>;