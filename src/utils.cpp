#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
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
	

void utils::print1DVector(const std::vector<unsigned>& a){
	for(unsigned i = 0; i < a.size(); i++){
			std::cout << a[i] << std::endl;
	}
	return;
}


double utils::vectors_distance(const std::vector<double>& a, const std::vector<double>& b)
{		
		std::vector<double> temp;
		
		if(b.size() == 0)
			return 0;
		
		std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp),[](double e1, double e2) {return pow((e1-e2),2);});
	
		return sqrt(std::accumulate(temp.begin(), temp.end(), 0.0));
}

// Find the intersect of two vectors
std::vector<unsigned> utils::intersect(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted){
	std::vector<unsigned> ret;
	
	if(!isSorted){
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}
	
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(ret));
	
	return ret;
}
