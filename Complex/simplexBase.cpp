#include <string>
#include <vector>
#include <typeinfo>
#include "simplexBase.hpp"
#include "simplexTree.hpp"
#include "simplexArrayList.hpp"
#include "indSimplexTree.hpp"

simplexBase::simplexBase(){return;}

simplexBase::simplexBase(std::map<std::string, std::string> configMap){
	setConfig(configMap);
	
	return;
}

simplexBase::simplexBase(double maxE, int maxDim){
	maxEpsilon = maxE;	
	maxDimension = maxDim;	
}

void simplexBase::setConfig(std::map<std::string, std::string> configMap){
	std::string debug;
	std::string outputFile;
	
	auto pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = std::atoi(configMap["outputFile"].c_str());
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end())
		maxDimension = std::atoi(configMap["dimensions"].c_str());
	else return;
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return;	
	
	std::cout << "Setting utils for : " << simplexType << std::endl;
	ut = utils(debug, outputFile);
	
	return;
}


void simplexBase::setDistanceMatrix(std::vector<std::vector<double>> _distMatrix){
	distMatrix = _distMatrix;
	return;
}

// simplexTree constructor, currently no needed information for the class constructor
simplexBase* simplexBase::newSimplex(const std::string &simplexT, std::map<std::string, std::string> configMap){
	simplexType = simplexT;
	
	if(simplexType == "simplexTree"){
		auto t = new simplexTree(maxEpsilon, distMatrix, maxDimension);
		t->setConfig(configMap);
		return t;
	} else if (simplexType == "simplexArrayList"){	
		auto t = new simplexArrayList(maxEpsilon, maxDimension, distMatrix);
		t->setConfig(configMap);
		return t;
	} else if (simplexType == "indSimplexTree"){
		auto t = new indSimplexTree(maxEpsilon, distMatrix, maxDimension);
		t->setConfig(configMap);
		return t;
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

void simplexBase::insertIterative(std::vector<double>&){
	ut.writeLog(simplexType,"No insert iterative function defined");
	return;
}

void simplexBase::deleteIterative(int){
	ut.writeLog(simplexType,"No delete iterative function defined");
	return;
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
	
