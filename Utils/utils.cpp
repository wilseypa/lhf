#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <numeric>
#include <iostream>
#include "utils.hpp"

// utils constructor, currently no needed information for the class constructor
utils::utils(){

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

void utils::print1DSet(const auto& a){
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

std::vector<double> utils::feature_distance(std::vector<double>* a, std::vector<double>* b){
	std::vector<double> ret;
	
	return ret;
}

double utils::vectors_distance(const double& a, const double& b){
		return pow((a-b),2);
}

double utils::vectors_distance(const std::vector<double>& a, const std::vector<double>& b){		
		std::vector<double> temp;
		
		if(b.size() == 0)
			return 0;
		
		std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp),[](double e1, double e2) {return pow((e1-e2),2);});
	
		
	
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
	
	for(auto iter = v1.begin(); iter!= v1.end(); iter++){
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
	std::cout << std::endl;
	
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
std::vector<unsigned> utils::symmetricDiff(std::set<unsigned> v1, std::set<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;
	std::vector<unsigned> retTemp;
	
	if(v1 == v2)
		return ret;
	
	//if(!isSorted){
	//	sort(v1.begin(), v1.end());
	//	sort(v2.begin(), v2.end());
//	}
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

