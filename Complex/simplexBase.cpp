#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <typeinfo>
#include "simplexBase.hpp"
#include "simplexTree.hpp"
#include "simplexArrayList.hpp"
#include "alphaComplex.hpp"
#include "witnessComplex.hpp"
#include "betaComplex.hpp"

template<typename nodeType>
simplexBase<nodeType>::simplexBase(){return;}

template<typename nodeType>
simplexBase<nodeType>::simplexBase(std::map<std::string, std::string> &configMap){
	setConfig(configMap);
	return;
}

template<typename nodeType>
simplexBase<nodeType>::simplexBase(double maxE, int maxDim){
	maxEpsilon = maxE;
	maxDimension = maxDim;
}

template<typename nodeType>
void simplexBase<nodeType>::setConfig(std::map<std::string, std::string> &configMap){
	std::string debug;
	std::string outputFile;

	auto pipe = configMap.find("debug");
	if(pipe != configMap.end())
		debug = std::atoi(configMap["debug"].c_str());
	pipe = configMap.find("outputFile");
	if(pipe != configMap.end())
		outputFile = configMap["outputFile"].c_str();
	pipe = configMap.find("dimensions");
	if(pipe != configMap.end())
		maxDimension = std::atoi(configMap["dimensions"].c_str());
	else return;
	pipe = configMap.find("epsilon");
	if(pipe != configMap.end())
		maxEpsilon = std::atof(configMap["epsilon"].c_str());
	else return;
	pipe = configMap.find("complexType");
	if(pipe != configMap.end())
		complexType = configMap["complexType"];
	
	pipe = configMap.find("simplicialComplex");
	if(pipe != configMap.end())
		simplicialComplex = configMap["simplicialComplex"];
	
	pipe = configMap.find("alphaFilterationValue");
	if(pipe != configMap.end())
		alphaFilterationValue = std::atof(configMap["alphaFilterationValue"].c_str());
	else return;
	ut = utils(debug, outputFile);
	ut.writeLog(simplexType,"Configured utils for : " + simplexType);

	return;
}


template<typename nodeType>
void simplexBase<nodeType>::setDistanceMatrix(std::vector<std::vector<double>>* _distMatrix){
	distMatrix = _distMatrix;
	return;
}

template<typename nodeType>
void simplexBase<nodeType>::setIncidenceMatrix(std::vector<std::vector<bool>>* _incidenceMatrix){
	incidenceMatrix = _incidenceMatrix;
	return;
}
// simplexTree constructor, currently no needed information for the class constructor
template<typename nodeType>
simplexBase<nodeType>* simplexBase<nodeType>::newSimplex(const std::string &simplexT, std::map<std::string, std::string> &configMap){
    std::cout << "Setting up " << simplexT << std::endl;
    
	if(simplexT == "simplexTree"){
		auto t = new simplexTree<nodeType>(0, 0);
		t->setConfig(configMap);
		return t;
	} else if (simplexT == "simplexArrayList"){
		auto t = new simplexArrayList<nodeType>(0, 0);
		t->setConfig(configMap);
		return t;
	} else if (simplexT == "alphaComplex"){
		auto t = new alphaComplex<nodeType>(0, 0);
		t->setConfig(configMap);
		return t;
	} else if (simplexT == "witnessComplex"){
		auto t = new witnessComplex<nodeType>(0, 0);
		t->setConfig(configMap);
		return t;
	} else if (simplexT == "betaComplex"){
		auto t = new betaComplex<nodeType>(0, 0);
		t->setConfig(configMap);
		return t;
	}
	return 0;
}

template<typename nodeType>
void simplexBase<nodeType>::setEnclosingRadius(double r){
	maxEpsilon = r;
}

template<typename nodeType>
std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> simplexBase<nodeType>::getDimEdges(int dim){
	if(dim >= simplexList.size()){
		ut.writeLog(simplexType,"Error: requested dimension beyond complex");
		std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> a;
		return a;
	}
	return simplexList[dim];
}

template<typename nodeType>
std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> simplexBase<nodeType>::getdelaunayDimEdges(int dim){
	ut.writeLog(simplexType,"No getdelunayDimEdges function defined");
	std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>> simplexList[dim];
	return simplexList[dim];
}

template<typename nodeType>
std::vector<std::set<std::shared_ptr<nodeType>, cmpByWeight<std::shared_ptr<nodeType>>>> simplexBase<nodeType>::getAllEdges(){
	return simplexList;
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexBase<nodeType>::getAllCofacets(const std::set<unsigned>& simplex){
	return getAllCofacets(simplex, 0, std::unordered_map<std::shared_ptr<nodeType>, std::shared_ptr<nodeType>>(), false);
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexBase<nodeType>::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<std::shared_ptr<nodeType>, std::shared_ptr<nodeType>>& pivotPairs, bool checkEmergent){
	ut.writeLog(simplexType,"No get cofacets function defined");
	return std::vector<std::shared_ptr<nodeType>>();
}

template<typename nodeType>
std::vector<nodeType*> simplexBase<nodeType>::getAllCofacets(std::shared_ptr<nodeType>, const std::unordered_map<long long, std::shared_ptr<nodeType>>&, bool){
	ut.writeLog(simplexType,"No get cofacets function defined");
	return std::vector<nodeType*>();
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexBase<nodeType>::getAllDelaunayCofacets(std::shared_ptr<nodeType>){
	ut.writeLog(simplexType,"No getdelaunay cofacets function defined");
	return std::vector<std::shared_ptr<nodeType>>();
}

template<typename nodeType>
std::vector<nodeType*> simplexBase<nodeType>::getAllDelaunayCofacets_basePointer(std::shared_ptr<nodeType>){
	ut.writeLog(simplexType,"No getdelaunay cofacets function defined");
	return std::vector<nodeType*>();
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexBase<nodeType>:: getAllDelaunayCofacets(std::shared_ptr<nodeType> simp, std::unordered_map<std::shared_ptr<nodeType>,std::shared_ptr<nodeType>> pivotPairs,bool emergent){
	ut.writeLog(simplexType,"No getdelaunay cofacets function defined");
	return std::vector<std::shared_ptr<nodeType>>();
}
template<typename nodeType>
std::vector<nodeType*> simplexBase<nodeType>::getAllCofacets(std::shared_ptr<nodeType>){
	ut.writeLog(simplexType,"No get cofacets function defined");
	return std::vector<nodeType*>();
}

template<typename nodeType>
std::vector<nodeType*> simplexBase<nodeType>::getAllFacets(nodeType*){
	ut.writeLog(simplexType,"No get facets function defined");
	return std::vector<nodeType*>();
}

template<typename nodeType>
std::vector<nodeType*> simplexBase<nodeType>::getAllFacets(std::shared_ptr<nodeType> simplex){
	return getAllFacets(simplex.get());
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexBase<nodeType>::getAllFacets_P(std::shared_ptr<nodeType> simplex){
	ut.writeLog(simplexType,"No get facets function defined");
	return std::vector<std::shared_ptr<nodeType>>();
}

template<typename nodeType>
double simplexBase<nodeType>::getSize(){
	ut.writeLog(simplexType,"No size function defined");
	return -1;
}

template<typename nodeType>
bool simplexBase<nodeType>::insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&){
	ut.writeLog(simplexType,"No insert iterative function defined");
	return false;
}

template<typename nodeType>
bool simplexBase<nodeType>::insertIterative(std::vector<double>&, std::vector<std::vector<double>>&){
	ut.writeLog(simplexType,"No insert iterative function defined");
	return false;
}


template<typename nodeType>
void simplexBase<nodeType>::deleteIterative(int){
	ut.writeLog(simplexType,"No delete iterative function defined");
	return;
}


template<typename nodeType>
void simplexBase<nodeType>::deleteIndexRecurse(int){
	ut.writeLog(simplexType,"No recursive delete function defined");
	return;
}


template<typename nodeType>
void simplexBase<nodeType>::insert(){
	ut.writeLog(simplexType,"No insert function defined");
	return;
}

template<typename nodeType>
bool simplexBase<nodeType>::find(std::vector<unsigned>){
	ut.writeLog(simplexType,"No find function defined");
	return false;
}

template<typename nodeType>
bool simplexBase<nodeType>::find(std::set<unsigned>){
	ut.writeLog(simplexType,"No find function defined");
	return false;
}

template<typename nodeType>
int simplexBase<nodeType>::vertexCount(){
	ut.writeLog(simplexType,"No vertexCount function defined");
	return -1;
}

template<typename nodeType>
void simplexBase<nodeType>::prepareCofacets(int dim){
	ut.writeLog(simplexType,"No prepareCofacets function defined");
	return;
}

template<typename nodeType>
void simplexBase<nodeType>::prepareFacets(int dim){
	ut.writeLog(simplexType,"No prepareFacets function defined");
	return;
}

template<typename nodeType>
int simplexBase<nodeType>::simplexCount(){
	ut.writeLog(simplexType,"No simplexCount function defined");
	return -1;
}

template<typename nodeType>
void simplexBase<nodeType>::outputComplex(){
	ut.writeLog(simplexType,"No outputComplex function defined");
	return;
}


template<typename nodeType>
void simplexBase<nodeType>::expandDimensions(int dim){
	ut.writeLog(simplexType,"No expandDimensions function defined");
	return;
}

/*template<typename nodeType>
void simplexBase<nodeType>::reduceComplex(){
	ut.writeLog(simplexType,"No reduceComplex function defined");
	return;
}*/

template<typename nodeType>
void simplexBase<nodeType>::setStreamEvaluator(bool (*f) (std::vector<double>&, std::vector<std::vector<double>>&)){
	streamEval = f;
	ut.writeLog(simplexType,"Changed stream evaluator");
	return;
}


template<typename nodeType>
bool simplexBase<nodeType>::streamEvaluator(std::vector<double>& vector, std::vector<std::vector<double>>& window){
	//Do some evaluation of whether the point should stay or not
	//		For now, let's look at the deviation of connections

	auto reps = ut.nearestNeighbors(vector, window);

	double sum = std::accumulate(reps.begin(), reps.end(), 0.0);
	double mean = sum / reps.size();

	std::vector<double> diff(reps.size());
	std::transform(reps.begin(), reps.end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / reps.size());

	stats += std::to_string(runningVectorCount) + "," + std::to_string(mean) + "," + std::to_string(stdev) + ",";

	std::sort(reps.begin(), reps.end());
	std::vector<double> kNN;
	int k = 20;

	for(int i = 0; i < k; i++){
		kNN.push_back(reps[i]);
	}

	double sum_NN = std::accumulate(kNN.begin(), kNN.end(), 0.0);
	double mean_NN = sum_NN / kNN.size();

	std::vector<double> diff_NN(kNN.size());
	std::transform(kNN.begin(), kNN.end(), diff_NN.begin(),std::bind2nd(std::minus<double>(), mean_NN));
	double sq_sum_NN = std::inner_product(diff_NN.begin(), diff_NN.end(), diff_NN.begin(), 0.0);
	double stdev_NN = std::sqrt(sq_sum_NN / kNN.size());

	stats += std::to_string(k) + "," + std::to_string(mean_NN) + "," + std::to_string(stdev_NN) + ",";

	if (true){//stdev_NN > 10000){
		//std::cout << "\tAccept: (stdev > 0.5 , " << stdev << ")" << std::endl;
		stats += "Accept\n";
		return true;
	}
	stats += "Reject\n";
	//std::cout << "\tReject: (stdev > 0.5 , " << stdev << ")" << std::endl;
	return false;
}

template<typename nodeType>
simplexBase<nodeType>::~simplexBase(){}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexBase<nodeType>::expandDimension(std::vector<std::shared_ptr<nodeType>> edges){
	std::vector<std::shared_ptr<nodeType>> ret;
	ut.writeLog(simplexType,"No expandDimension function defined");
	return ret;
}

template<typename nodeType>
std::vector<std::shared_ptr<nodeType>> simplexBase<nodeType>::expanddelaunayDimension(int dim){
	std::vector<std::shared_ptr<nodeType>> ret;
	ut.writeLog(simplexType,"No expanddelaunayDimension function defined");
	return ret;
}

template class simplexBase<simplexNode>;
template class simplexBase<alphaNode>;
template class simplexBase<witnessNode>;
