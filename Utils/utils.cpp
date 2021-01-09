#include <string>
#include <cmath>
#include <algorithm>
#include <utility>
#include <numeric>
#include <iostream>
#include <fstream>
#include "utils.hpp"



double d_simplexIntersection:: Det2D(TriPoint &p1, TriPoint &p2, TriPoint &p3)
{
	return +p1.first*(p2.second-p3.second)
		+p2.first*(p3.second-p1.second)
		+p3.first*(p1.second-p2.second);
}

void d_simplexIntersection::CheckTriWinding(TriPoint &p1, TriPoint &p2, TriPoint &p3, bool allowReversed)
{
	double detTri = Det2D(p1, p2, p3);
	if(detTri < 0.0)
	{
		if (allowReversed)
		{
			TriPoint a = p3;
			p3 = p2;
			p2 = a;
		}
		else throw std::runtime_error("triangle has wrong winding direction");
	}
}

bool d_simplexIntersection::BoundaryCollideChk(TriPoint &p1, TriPoint &p2, TriPoint &p3, double eps)
{
	return Det2D(p1, p2, p3) < eps;
}

bool d_simplexIntersection::BoundaryDoesntCollideChk(TriPoint &p1, TriPoint &p2, TriPoint &p3, double eps)
{
	return Det2D(p1, p2, p3) <= eps;
}

bool d_simplexIntersection::TriTri2D(TriPoint *t1,	TriPoint *t2,double eps = 0.0, bool allowReversed = true, bool onBoundary = false)
{
	//Trangles must be expressed anti-clockwise
	CheckTriWinding(t1[0], t1[1], t1[2], allowReversed);
	CheckTriWinding(t2[0], t2[1], t2[2], allowReversed);

	bool (*chkEdge)(TriPoint &, TriPoint &, TriPoint &, double) = NULL;
	if(onBoundary) //Points on the boundary are considered as colliding
		chkEdge = BoundaryCollideChk;
	else //Points on the boundary are not considered as colliding
		chkEdge = BoundaryDoesntCollideChk;

	//For edge E of trangle 1,
	for(unsigned i=0; i<3; i++)
	{
		unsigned j=(i+1)%3;

		//Check all Points of trangle 2 lay on the external side of the edge E. If
		//they do, the triangles do not collide.
		if (chkEdge(t1[i], t1[j], t2[0], eps) &&
			chkEdge(t1[i], t1[j], t2[1], eps) &&
			chkEdge(t1[i], t1[j], t2[2], eps))
			return false;
	}

	//For edge E of trangle 2,
	for(unsigned i=0; i<3; i++)
	{
		unsigned j=(i+1)%3;

		//Check all Points of trangle 1 lay on the external side of the edge E. If
		//they do, the triangles do not collide.
		if (chkEdge(t2[i], t2[j], t1[0], eps) &&
			chkEdge(t2[i], t2[j], t1[1], eps) &&
			chkEdge(t2[i], t2[j], t1[2], eps))
			return false;
	}

	//The triangles collide
	return true;
}




unionFind::unionFind(int n) : rank(n, 0), parent(n, 0) {
	for(int i=0; i<n; i++) parent[i]=i;
}

int unionFind::find(int i){
	if(i == parent[i]) return i; //Found name of the component
	parent[i] = find(parent[i]); //Path Compression
	return parent[i];
}

bool unionFind::join(int x, int y){ //Union by rank
	x = find(x);
	y = find(y);
	if(x == y) return false;
	if(rank[x] == rank[y]){
		rank[y]++;
		parent[x] = y;
	} else if(rank[x] < rank[y]){
		parent[x] = y;
	} else{
		parent[y] = x;
	}
	return true;
}


// utils constructor, currently no needed information for the class constructor
utils::utils(){}

utils::utils(std::string _debug, std::string _outputFile){
	debug = _debug;
	outputFile = _outputFile;
}

double utils::computeMaxRadius(int k, std::vector<std::vector<double>> &centroids, std::vector<std::vector<double>> &originalData, std::vector<unsigned> &labels){
	double maxRadius = 0;
	double curRadius = 0;
	
	//Iterate through each point
	for(unsigned i = 0; i < originalData.size(); i++){
		
		//Check the distance of this point to it's centroid
		curRadius = vectors_distance(originalData[i], centroids[labels[i]]);
		
		if(curRadius > maxRadius)
			maxRadius = curRadius;
	}
	
	return maxRadius;	
}

double utils::computeAvgRadius(int k, std::vector<std::vector<double>> &centroids, std::vector<std::vector<double>> &originalData, std::vector<unsigned> &labels){
	double totalRadius = 0;
	
	//Iterate through each point
	for(unsigned i = 0; i < originalData.size(); i++){
		
		//Check the distance of this point to it's centroid
		totalRadius += vectors_distance(originalData[i], centroids[labels[i]]);
	}
	
	return totalRadius / originalData.size();	
	
}

std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> utils::separatePartitions(int k, std::vector<std::vector<double>> originalData, std::vector<unsigned> labels){
	std::vector<std::vector<double>> a;
	std::vector<unsigned> b;
	std::vector<std::vector<std::vector<double>>> res(k, a);
	std::vector<std::vector<unsigned>> labres(k, b);
	
	for(unsigned i = 0; i < labels.size(); i++){
		res[labels[i]].push_back(originalData[i]);
		labres[labels[i]].push_back(i);
	}
	
	return std::make_pair(labres, res);
}

std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> utils::separatePartitions(double rad, std::vector<std::vector<double>> centroids, std::vector<std::vector<double>> originalData, std::vector<unsigned> labels){
	std::vector<std::vector<double>> a;
	std::vector<unsigned> b;
	
	//Store the results to return
	std::vector<std::vector<std::vector<double>>> res(centroids.size(), a);
	std::vector<std::vector<unsigned>> labres(centroids.size(), b);

	//Iterate through each label
	for(unsigned i = 0; i < labels.size(); i++){
		
		//Check for this point belonging to each centroid
		for(unsigned j = 0; j < centroids.size(); j++){
			if(labels[i] == j){
				//If this is a labeled constituent put to front of array
				res[j].insert(res[j].begin(), originalData[i]);
				labres[j].insert(labres[j].begin(), i);
				
			}else{
				
				double curRad = vectors_distance(originalData[i], centroids[j]);
				
				//If this distance is less than our radius cutoff push to back
				if(curRad < rad){
					res[j].push_back(originalData[i]);
					labres[j].push_back(i);	
				}	
			}
		}
	}
	
	return std::make_pair(labres, res);
}

std::vector<std::vector<std::vector<double>>> utils::separateBoundaryPartitions(std::vector<std::set<unsigned>> boundaryLists, std::vector<std::vector<double>> originalData, std::vector<unsigned> labels){
	std::vector<std::vector<double>> a;
	std::vector<std::vector<std::vector<double>>> res(boundaryLists.size(), a);
		
	for(unsigned i = 0; i < originalData.size(); i++){
		
		for(unsigned j = 0; j < boundaryLists.size(); j++){
			
			if(boundaryLists[j].find(labels[i]) != boundaryLists[j].end())
				res[j].push_back(originalData[i]);
		}
	}
	
	return res;
}

// void utils::extractBoundaryPoints(std::vector<bettiBoundaryTableEntry>& bettiTable){
// 	for(auto& bet: bettiTable){
// 		std::set<unsigned> bound;
// 		for(auto simplex : bet.boundary) bound.insert(simplex->simplex.begin(), simplex->simplex.end());
// 		bet.boundaryPoints = bound;
// 	}
// }

std::set<unsigned> utils::extractBoundaryPoints(std::vector<simplexNode_P> boundary){
	std::set<unsigned> boundaryPoints;
	for(auto simplex : boundary) boundaryPoints.insert(simplex->simplex.begin(), simplex->simplex.end());
	return boundaryPoints;
}

std::vector<bettiBoundaryTableEntry> utils::mapPartitionIndexing(std::vector<unsigned> partitionedLabels, std::vector<bettiBoundaryTableEntry> bettiTable){
	for(auto& bet : bettiTable){
		std::set<unsigned> convBound;
		
		for(auto ind : bet.boundaryPoints){
			convBound.insert(partitionedLabels[ind]);
		}
		
		bet.boundaryPoints = convBound;
	}
	return bettiTable;
}

void utils::print2DVector(const std::vector<std::vector<unsigned>>& a){
	for(unsigned i = 0; i < a.size(); i++){
		for(unsigned j = 0; j < a[i].size(); j++){
			std::cout << a[i][j] << '\t';
		}
		std::cout << std::endl;
	}
	return;
}

void utils::print1DSet(const std::pair<std::set<unsigned>, double>& a){
		std::cout << "Test\t";
		
		for(auto iter = a.first.begin(); iter!= a.first.end(); iter++){
			std::cout << *iter << ",";
		}
		std::cout << "\t";
}

void utils::print1DVector(const std::vector<double>& a){
	for(unsigned i = 0; i < a.size(); i++){
			std::cout << a[i] << ",";
	}
	std::cout << "\n";
	return;
}

void utils::print1DVector(const std::vector<unsigned>& a){
	for(unsigned i = 0; i < a.size(); i++){
			std::cout << a[i] << ",";
	}
	std::cout << "\n";
	return;
}

void utils::print1DVector(const std::set<unsigned>& a){
	for(auto z : a){
			std::cout << z << ",";
	}
	std::cout << "\n";
	return;
}

double utils::vectors_distance(const double& a, const double& b){
		return pow((a-b),2);
}

std::set<unsigned> utils::setXOR(std::set<unsigned>& setA, std::set<unsigned>& setB){
	std::set<unsigned> ret;

	//if(setA.size() == 0)
	//	ret = setB;
	//else if (setB.size() == 0)
	//	ret = setA;
	//else
	set_symmetric_difference(setA.begin(), setA.end(), setB.begin(), setB.end(), std::inserter(ret, ret.begin()));
	
	return ret;
}

double utils::vectors_distance(const std::vector<double>& a, const std::vector<double>& b){		
	std::vector<double> temp;
	
	if(b.size() == 0)
		return 0;
	
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp),[](double e1, double e2) {
		return pow((e1-e2),2);
	});

	return sqrt(std::accumulate(temp.begin(), temp.end(), 0.0));
}

std::vector<unsigned> utils::setIntersect(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;
	
	if(v1 == v2)
		return v1;
	
	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}
	
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(ret));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = ret.begin(); iter!= ret.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/
	
	return ret;
	
}

std::set<unsigned> utils::setIntersect(std::set<unsigned> v1, std::set<unsigned> v2, bool isSorted){
	std::set<unsigned> ret;
	
	if(v1 == v2)
		return v1;
	
	//if(!isSorted){
	//	sort(v1.begin(), v1.end());
	//	sort(v2.begin(), v2.end());
	//}
	
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(ret, ret.begin()));
	
	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = ret.begin(); iter!= ret.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/
	
	return ret;
	
}
	


// Find the intersect of two vectors
std::pair<std::vector<unsigned>, std::vector<unsigned>> utils::intersect(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::pair<std::vector<unsigned>, std::vector<unsigned>> ret;
	std::pair<std::vector<unsigned>, std::vector<unsigned>> retTemp;
	
	if(v1 == v2)
		return ret;
	
	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}
	
	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp.second));
	set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp.first));
	
	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.first.begin(); iter!= retTemp.first.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.second.begin(); iter!= retTemp.second.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;
	std::cout << std::to_string(retTemp.first.size() == 2) << std::endl;

	*/
	
	if(retTemp.first.size() == 2)
		return retTemp;
	else
		return ret;
}

// Find the symmetric difference of two vectors
std::vector<unsigned> utils::symmetricDiff(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;
	std::vector<unsigned> retTemp;
	
	if(v1 == v2)
		return ret;
	
	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}
	set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp));
	
	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.begin(); iter!= retTemp.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/
	
	return retTemp;
}

// Find the symmetric difference of two vectors
std::set<unsigned> utils::symmetricDiff(std::set<unsigned> v1, std::set<unsigned> v2, bool isSorted){
	std::set<unsigned> ret;
	std::set<unsigned> retTemp;
	
	if(v1 == v2)
		return ret;
	
	//if(!isSorted){
	//	sort(v1.begin(), v1.end());
	//	sort(v2.begin(), v2.end());
//	}
	set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(),  std::inserter(retTemp, retTemp.begin()));
	
	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.begin(); iter!= retTemp.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/
	
	return retTemp;
}

//Iteratively build subsets (faces) of the simplex set
std::vector<std::set<unsigned>> utils::getSubsets(std::set<unsigned> set, int dim){
	std::vector<std::set<unsigned>> subset;
	std::set<unsigned> empty;
	subset.push_back(empty);

	//For each set in the 
	for(auto i = set.begin(); i!= set.end(); i++){
		std::vector<std::set<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++){
			subsetTemp[j].insert(entry);
		}
		
		unsigned z = 0;
		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++){
			subset.push_back(*j);
			
		}
	}
	
	std::vector<std::set<unsigned>> retSubset;
	
	for(std::set<unsigned> z : subset){
		if(z.size() == dim)
			retSubset.push_back(z);
	}
	return retSubset;
}

// Find the union of two vectors
std::set<unsigned> utils::setUnion(std::set<unsigned> v1, std::set<unsigned> v2){
	std::set<unsigned> retTemp;
	
	if(v1 == v2) return v1;
	
	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(retTemp, retTemp.begin()));
		
	return retTemp;
}



// Find the union of two vectors
std::vector<unsigned> utils::setUnion(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;
	std::vector<unsigned> retTemp;
	
	if(v1 == v2)
		return ret;
	
	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}
	
	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp));
		
	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.begin(); iter!= retTemp.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/
	return retTemp;
}

void utils::writeLog(std::string module, std::string message){
	if(debug == "1" || debug == "true"){
		std::cout << "[" << module << "]:\t" << message << std::endl;
	} else {
		writeFile("[" + module + "]:\t" + message);
	}
	return;
}

void utils::writeDebug(std::string module, std::string message){
	if(debug == "0" || debug == "false"){
		return;
	} else {
		std::cout << "[DEBUG]\t[" << module << "]:\t" << message << std::endl;
		writeFile("[DEBUG]\t[" + module + "]:\t" + message);
	}
	
	return;
}

void utils::writeFile(std::string fullMessage){
	std::ofstream outfile;
	outfile.open(outputFile+"_debug.txt", std::ios_base::app);
	outfile << fullMessage << "\n"; 
	
	return;
}


bool utils::sortBySecond(const std::pair<std::set<unsigned>, double> &a, const std::pair<std::set<unsigned>, double> &b){
	if(a.second == b.second){ //If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
		auto itA = a.first.rbegin(), itB = b.first.rbegin();
		while(itA != a.first.rend()){
			if(*itA != *itB) return *itA>*itB;
			++itA; ++itB;
		}
		return false;
	} else{
		return a.second<b.second;
	}
}

//Iteratively build subsets (faces) of the simplex set
std::vector<std::set<unsigned>> utils::getSubsets(std::set<unsigned> set){
	std::vector<std::set<unsigned>> subset;
	std::set<unsigned> empty;
	subset.push_back(empty);

	//For each set in the 
	for(auto i = set.begin(); i!= set.end(); i++){
		std::vector<std::set<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++){
			subsetTemp[j].insert(entry);
		}
		
		unsigned z = 0;
		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++){
			subset.push_back(*j);
			
		}
	}
	
	std::vector<std::set<unsigned>> retSubset;
	
	for(std::set<unsigned> z : subset){
		if(z.size() == set.size() - 1)
			retSubset.push_back(z);
	}
	return retSubset;
}

//Iteratively build subsets (faces) of the simplex set
std::vector<std::vector<unsigned>> utils::getSubsets(std::vector<unsigned> set){
	std::vector<std::vector<unsigned>> subset;
	std::vector<unsigned> empty;
	subset.push_back(empty);

	//For each set in the 
	for(auto i = set.begin(); i!= set.end(); i++){
		std::vector<std::vector<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++){
			subsetTemp[j].push_back(entry);
		}
		
		unsigned z = 0;
		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++){
			subset.push_back(*j);
			
		}
	}
	
	std::vector<std::vector<unsigned>> retSubset;
	
	for(std::vector<unsigned> z : subset){
		if(z.size() == set.size() - 1)
			retSubset.push_back(z);
	}
	return retSubset;
}


std::vector<double> utils::nearestNeighbors(std::vector<double>& point, std::vector<std::vector<double>>& pointcloud){
	//based on random projection, x is current point being examined, n is number of centroids/facilities
	utils ut;
	std::vector<double> retVal;
		
	//Get sq distances for each point
	for(auto currentPoint : pointcloud){
		retVal.push_back(ut.vectors_distance(point, currentPoint));
	}
	
	return retVal;

}


std::vector<std::vector<double>> utils::deserialize(std::vector<double> serialData, unsigned dim){
	
	//First check if the vector size matches the dimension
	if(serialData.size() % dim != 0){
		std::cout << "Error occurred when deserializing data: invalid size" << std::endl;
		return {};
	}
	
	//Deduce the number of vectors
	auto n = serialData.size() / dim;
	auto begin = serialData.begin();
	
	std::vector<std::vector<double>> ret(n, std::vector<double>(dim));
	
	for(unsigned i = 0; i < n; i++){
		for(unsigned j = 0; j < dim; j++){
			ret[i][j] = (serialData)[(i*dim) + j];
		}
	}
	return ret;
}


std::vector<double> utils::serialize(std::vector<std::vector<double>> &origData){
	
	//Make sure we have data to serialize
	if(origData.size() == 0){
		std::cout << "Error occurred when serializing data: empty vector argument" << std::endl;
		return {};
	}
	
	auto n = origData.size();
	auto d = origData[0].size();
	
	//Preallocate our vector to prevent resizing
	std::vector<double> ret(n*d);
	
	//Copy element by element
	for(unsigned i = 0; i < n; i++){
		for(unsigned k = 0; k < d; k++){
			ret[(i*d)+k] = origData[i][k];
		}
	}
	
	return ret;
}
