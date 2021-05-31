#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <typeinfo>
#include "simplexBase.hpp"
//#include "simplexTree.hpp"
#include "simplexArrayList.hpp"

template<typename T>
simplexBase<T>::simplexBase(){return;}

template<typename T>
simplexBase<T>::simplexBase(std::map<std::string, std::string> &configMap){
	setConfig(configMap);
	return;
}

template<typename T>
simplexBase<T>::simplexBase(double maxE, int maxDim){
	maxEpsilon = maxE;
	maxDimension = maxDim;
}

template<typename T>
void simplexBase<T>::setConfig(std::map<std::string, std::string> &configMap){
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
		simplexType = configMap["complexType"];
	
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


template<typename T>
void simplexBase<T>::setDistanceMatrix(std::vector<std::vector<double>>* _distMatrix){
	distMatrix = _distMatrix;
	return;
}

// simplexTree constructor, currently no needed information for the class constructor
template<typename T>
simplexBase<T>* simplexBase<T>::newSimplex(const std::string &simplexT, std::map<std::string, std::string> &configMap){
	if(simplexT == "simplexTree"){
		//maxEpsilon and maxDimension are overwritten by setConfig
		//auto t = new simplexTree(0, 0);
		//t->setConfig(configMap);
		auto t = new simplexArrayList<T>(0, 0);
		t->setConfig(configMap);
		return t;
	} else if (simplexT == "simplexArrayList"){
		auto t = new simplexArrayList<T>(0, 0);
		t->setConfig(configMap);
		return t;
	}
	return 0;
}

template<typename T>
void simplexBase<T>::setEnclosingRadius(double r){
	maxEpsilon = r;
}

template<typename T>
std::set<std::shared_ptr<T>, cmpByWeight<std::shared_ptr<T>>> simplexBase<T>::getDimEdges(int dim){
	if(dim >= simplexList.size()){
		ut.writeLog(simplexType,"Error: requested dimension beyond complex");
		std::set<std::shared_ptr<T>, cmpByWeight<std::shared_ptr<T>>> a;
		return a;
	}
	return simplexList[dim];
}

template<typename T>
std::vector<std::set<std::shared_ptr<T>, cmpByWeight<std::shared_ptr<T>>>> simplexBase<T>::getAllEdges(){
	return simplexList;
}

template<typename T>
std::vector<std::shared_ptr<T>> simplexBase<T>::getAllCofacets(const std::set<unsigned>& simplex){
	return getAllCofacets(simplex, 0, std::unordered_map<std::shared_ptr<T>, std::shared_ptr<T>>(), false);
}

template<typename T>
std::vector<std::shared_ptr<T>> simplexBase<T>::getAllCofacets(const std::set<unsigned>& simplex, double simplexWeight, const std::unordered_map<std::shared_ptr<T>, std::shared_ptr<T>>& pivotPairs, bool checkEmergent){
	ut.writeLog(simplexType,"No get cofacets function defined");
	return std::vector<std::shared_ptr<T>>();
}

template<typename T>
std::vector<T*> simplexBase<T>::getAllCofacets(std::shared_ptr<T>, const std::unordered_map<long long, std::shared_ptr<T>>&, bool){
	ut.writeLog(simplexType,"No get cofacets function defined");
	return std::vector<T*>();
}

template<typename T>
std::vector<std::shared_ptr<T>> simplexBase<T>::getAllDelaunayCofacets(std::shared_ptr<T>){
	ut.writeLog(simplexType,"No getdelaunay cofacets function defined");
	return std::vector<std::shared_ptr<T>>();
}

template<typename T>
std::vector<T*> simplexBase<T>::getAllCofacets(std::shared_ptr<T>){
	ut.writeLog(simplexType,"No get cofacets function defined");
	return std::vector<T*>();
}

template<typename T>
std::vector<T*> simplexBase<T>::getAllFacets(T*){
	ut.writeLog(simplexType,"No get facets function defined");
	return std::vector<T*>();
}

template<typename T>
std::vector<T*> simplexBase<T>::getAllFacets(std::shared_ptr<T> simplex){
	return getAllFacets(simplex.get());
}

template<typename T>
std::vector<std::shared_ptr<T>> simplexBase<T>::getAllFacets_P(std::shared_ptr<T> simplex){
	ut.writeLog(simplexType,"No get facets function defined");
	return std::vector<std::shared_ptr<T>>();
}

template<typename T>
double simplexBase<T>::getSize(){
	ut.writeLog(simplexType,"No size function defined");
	return -1;
}

template<typename T>
bool simplexBase<T>::insertIterative(std::vector<double>&, std::vector<std::vector<double>>&, int&, int&){
	ut.writeLog(simplexType,"No insert iterative function defined");
	return false;
}

template<typename T>
bool simplexBase<T>::insertIterative(std::vector<double>&, std::vector<std::vector<double>>&){
	ut.writeLog(simplexType,"No insert iterative function defined");
	return false;
}


template<typename T>
void simplexBase<T>::deleteIterative(int){
	ut.writeLog(simplexType,"No delete iterative function defined");
	return;
}


template<typename T>
void simplexBase<T>::deleteIndexRecurse(int){
	ut.writeLog(simplexType,"No recursive delete function defined");
	return;
}


template<typename T>
void simplexBase<T>::insert(){
	ut.writeLog(simplexType,"No insert function defined");
	return;
}

template<typename T>
bool simplexBase<T>::find(std::vector<unsigned>){
	ut.writeLog(simplexType,"No find function defined");
	return false;
}

template<typename T>
bool simplexBase<T>::find(std::set<unsigned>){
	ut.writeLog(simplexType,"No find function defined");
	return false;
}

template<typename T>
int simplexBase<T>::vertexCount(){
	ut.writeLog(simplexType,"No vertexCount function defined");
	return -1;
}

template<typename T>
void simplexBase<T>::prepareCofacets(int dim){
	ut.writeLog(simplexType,"No prepareCofacets function defined");
	return;
}

template<typename T>
void simplexBase<T>::prepareFacets(int dim){
	ut.writeLog(simplexType,"No prepareFacets function defined");
	return;
}

template<typename T>
int simplexBase<T>::simplexCount(){
	ut.writeLog(simplexType,"No simplexCount function defined");
	return -1;
}

template<typename T>
void simplexBase<T>::outputComplex(){
	ut.writeLog(simplexType,"No outputComplex function defined");
	return;
}


template<typename T>
void simplexBase<T>::expandDimensions(int dim){
	ut.writeLog(simplexType,"No expandDimensions function defined");
	return;
}

/*template<typename T>
void simplexBase<T>::reduceComplex(){
	ut.writeLog(simplexType,"No reduceComplex function defined");
	return;
}*/

template<typename T>
void simplexBase<T>::setStreamEvaluator(bool (*f) (std::vector<double>&, std::vector<std::vector<double>>&)){
	streamEval = f;
	ut.writeLog(simplexType,"Changed stream evaluator");
	return;
}


template<typename T>
bool simplexBase<T>::streamEvaluator(std::vector<double>& vector, std::vector<std::vector<double>>& window){
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

template<typename T>
simplexBase<T>::~simplexBase(){}

template<typename T>
std::vector<std::shared_ptr<T>> simplexBase<T>::expandDimension(std::vector<std::shared_ptr<T>> edges){
	std::vector<std::shared_ptr<T>> ret;
	ut.writeLog(simplexType,"No expandDimension function defined");
	return ret;
}

//Explicit Template Class Instantiation
template class simplexBase<simplexNode>;
template class simplexBase<alphaNode>;