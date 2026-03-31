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
  //Yogesh will write evolutionary algorithm here
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

template class directPolytopal<simplexNode>;
template class directPolytopal<alphaNode>;
template class directPolytopal<witnessNode>;
