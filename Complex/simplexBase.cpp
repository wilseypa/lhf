#include <string>
#include <vector>
#include <typeinfo>
#include "simplexBase.hpp"
#include "simplexTree.hpp"
#include "simplexArrayList.hpp"
#include "indSimplexTree.hpp"

simplexBase::simplexBase(){return;}
simplexBase::simplexBase(double maxE, int maxDim){
	maxEpsilon = maxE;	
	maxDimension = maxDim;	
}
simplexBase::simplexBase(std::map<std::string, std::string> configMap){
	std::string debug;
	std::string outputFile;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = std::atoi(configMap["outputFile"].c_str());
	
	ut = utils(debug, outputFile);
	
	return;
}

void simplexBase::setDistanceMatrix(std::vector<std::vector<double>> _distMatrix){
	distMatrix = _distMatrix;
	return;
}

// simplexTree constructor, currently no needed information for the class constructor
simplexBase* simplexBase::newSimplex(const std::string &simplexT){
	simplexType = simplexT;
	
	if(simplexType == "simplexTree"){
		return new simplexTree(maxEpsilon, distMatrix, maxDimension);
	} else if (simplexType == "simplexArrayList"){	
		return new simplexArrayList(maxEpsilon, distMatrix);
	} else if (simplexType == "indSimplexTree"){
		return new indSimplexTree(maxEpsilon, distMatrix, maxDimension);
	}
	return 0;
}


std::vector<std::vector<unsigned>> simplexBase::getDimEdges(int dim, double epsilon){
	ut.writeLog("simplexBase","No get edges function defined for: " + simplexType);
	std::vector<std::vector<unsigned>> a;
	return a;
}

std::vector<std::vector<std::pair<std::set<unsigned>,double>>> simplexBase::getAllEdges(double epsilon){
	ut.writeLog("simplexBase","No get edges function defined for: " + simplexType);
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> a;
	return a;
}

std::vector<std::vector<indSimplexTree::graphEntry>> simplexBase::getIndexEdges(double epsilon){
	ut.writeLog("simplexBase","No get index edges function defined for: " + simplexType);
	std::vector<std::vector<indSimplexTree::graphEntry>> a;
	return a;
}

double simplexBase::getSize(){
	ut.writeLog("simplexBase","No size function defined for: " + simplexType);
	return 0;
}

void simplexBase::insert(std::vector<double>&){
	ut.writeLog("simplexBase","No insert function defined for: " + simplexType);
	return;
}

void simplexBase::find(std::vector<double>){
	ut.writeLog("simplexBase","No find function defined for: " + simplexType);
	return;
}

unsigned simplexBase::find(std::set<unsigned>){
	ut.writeLog("simplexBase","No find function defined for: " + simplexType);
	return -1;
}

int simplexBase::vertexCount(){
	ut.writeLog("simplexBase","No vertexCount function defined for: " + simplexType);
	return -1;
}

int simplexBase::simplexCount(){
	ut.writeLog("simplexBase","No simplexCount function defined for: " + simplexType);
	return -1;
}

void simplexBase::outputSimplex(){
	ut.writeLog("simplexBase","No outputSimplex function defined for: " + simplexType);
	return;
}

void simplexBase::expandDimensions(int dim){
	ut.writeLog("simplexBase","No expandDimensions function defined for: " + simplexType);
	return;
}
	
void simplexBase::reduceComplex(){
	ut.writeLog("simplexBase","No reduceComplex function defined for: " + simplexType);
	return;
}
	
