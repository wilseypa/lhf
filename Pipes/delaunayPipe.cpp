#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include "delaunayPipe.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include <fstream>

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

  int temp_id=0;
  int dim = inData.inputData[0].size();
  int no_of_points=inData.inputData.size();
  typedef CGAL::Delaunay_triangulation<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>     T;
  T dt(dim);
  auto point=T::Point (inData.inputData[0].begin(), inData.inputData[0].end());
  T::Vertex_handle hint;
  std::map<std::string, int> id;
  std::stringstream buffer;
  hint = dt.insert(point);
  buffer<<point;
  id.insert(std::make_pair(buffer.str(),temp_id++));
    for (int i=1;i<no_of_points;i++) {
      point=T::Point (inData.inputData[i].begin(), inData.inputData[i].end());
      hint = dt.insert(point, hint);
      std::stringstream().swap(buffer);
      buffer<<point;
      id.insert({buffer.str(),temp_id++});
  }
  std::vector<std::vector<unsigned>> dsimplexes;
  std::vector<unsigned> temp;
    for (auto i = dt.finite_full_cells_begin(); i != dt.finite_full_cells_end(); i++)
    {
      for (auto j = i->vertices_begin(); j != i->vertices_end(); j++)
      {
        std::stringstream().swap(buffer);
        buffer<<(**j);
        temp.push_back(id[buffer.str()]);
         // std::cout<< id[buffer.str()]<<" ";
      }
      dsimplexes.push_back(temp);
      temp.clear();
      //std::cout<<std::endl;
    }
  dt.clear();
  id.clear();
 ((alphaComplex<nodeType>*)inData.complex)->buildAlphaComplex(dsimplexes,inData.inputData.size(),inData.inputData);
  //std::cout<<dsimplexes.size();
  
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