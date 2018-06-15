#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "utils.hpp"

// utils constructor, currently no needed information for the class constructor
utils::utils(){

}


double utils::vectors_distance(const std::vector<double>& a, const std::vector<double>& b)
{		
		std::vector<double> temp;
		
		if(b.size() == 0)
			return 0;
		
		std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp),[](double e1, double e2) {return pow((e1-e2),2);});
	
		return sqrt(std::accumulate(temp.begin(), temp.end(), 0.0));
}

