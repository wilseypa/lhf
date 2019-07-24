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
	ut.writeLog(simplexType,"No get edges function defined");
	std::vector<std::vector<unsigned>> a;
	return a;
}

std::vector<std::vector<std::pair<std::set<unsigned>,double>>> simplexBase::getAllEdges(double epsilon){
	ut.writeLog(simplexType,"No get edges function defined");
	std::vector<std::vector<std::pair<std::set<unsigned>,double>>> a;
	return a;
}

std::vector<std::vector<indSimplexTree::graphEntry>> simplexBase::getIndexEdges(double epsilon){
	ut.writeLog(simplexType,"No get index edges function defined");
	std::vector<std::vector<indSimplexTree::graphEntry>> a;
	return a;
}

double simplexBase::getSize(){
	ut.writeLog(simplexType,"No size function defined");
	return 0;
}

void simplexBase::insert(std::vector<double>&){
	ut.writeLog(simplexType,"No insert function defined");
	return;
}

bool simplexBase::find(std::vector<unsigned>){
	ut.writeLog(simplexType,"No find function defined");
	return -1;
}

bool simplexBase::find(std::set<unsigned>){
	ut.writeLog(simplexType,"No find function defined");
	return -1;
}

int simplexBase::vertexCount(){
	ut.writeLog(simplexType,"No vertexCount function defined");
	return -1;
}

int simplexBase::simplexCount(){
	ut.writeLog(simplexType,"No simplexCount function defined");
	return -1;
}

void simplexBase::outputSimplex(){
	ut.writeLog(simplexType,"No outputSimplex function defined");
	return;
}

void simplexBase::expandDimensions(int dim){
	ut.writeLog(simplexType,"No expandDimensions function defined");
	return;
}
	
void simplexBase::reduceComplex(){
	ut.writeLog(simplexType,"No reduceComplex function defined");
	return;
}
	
