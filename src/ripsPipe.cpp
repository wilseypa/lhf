/*
 * ripsPipe hpp + cpp extend the basePipe class for calculating the 
 * VR complex from a neighborhood graph
 * 
 */

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <algorithm>
#include <set>
#include "ripsPipe.hpp"

// basePipe constructor
ripsPipe::ripsPipe(){
	pipeType = "Vietoris Rips (inductive)";
	return;
}


std::vector<unsigned> intersect(std::vector<unsigned> v1, std::vector<unsigned> v2){
	std::vector<unsigned> ret;
	
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(ret));
	
	return ret;
}


// runPipe -> Run the configured functions of this pipeline segment
pipePacket ripsPipe::runPipe(pipePacket inData){
	
	//Store the complex built from the VR expansion
	std::vector<std::vector<unsigned>> edges;
	std::vector<double> weights;
	std::vector<std::vector<unsigned>> temp;
	
	for(auto sortTemp : inData.workData.edges){
		sort(sortTemp.begin(), sortTemp.end());
		temp.push_back(sortTemp);
	}
	 
	//Iterate up to max dimension of simplex
	for(unsigned i = 1; i < dim; i++){
		std::vector<std::vector<unsigned>> test;
		
		
		//Iterate through each element in the previous dimension's edges
		for(unsigned j = 0; j < temp.size(); j++){
			std::vector<std::vector<unsigned>> tempIntersections;
			
			//First search for intersections of the current element
			for(unsigned t = j+1; t < temp.size(); t++){
				auto simp = intersect(temp[j], temp[t]);

				//This point intersects; potential candidate for a higher-level simplice
				if (simp.size() > 0){
					tempIntersections.push_back(temp[t]);
					
					//Loop through the current tempIntersections list
					for(unsigned z = 0; z < tempIntersections.size(); z++){
						simp = intersect(tempIntersections[z], temp[t]);
						
						if(simp.size() > 0){
							//Copy the tested simplex into the new intersection
							for(auto e : temp[j])
								simp.push_back(e);
								
							//remove any duplicates
							std::set<unsigned> s(simp.begin(), simp.end());
							simp.assign(s.begin(), s.end());
					
							test.push_back(simp);
						}
					
					
					}
				}
			}		
			
			/*//Search for intersections of elements that already intersected the current
			for(unsigned t = 0; t < tempIntersections.size(); t++){
				for(unsigned z = 0; z < t; z++){
					auto simp = intersect(tempIntersections[t],tempIntersections[z]);
	
					for(auto e : tempSimplex)
						simp.push_back(e);		
					
					
				}
			}	*/
		}		
		
		//Filter any duplicate simplices
		std::set<std::vector<unsigned>> s(test.begin(), test.end());
		test.assign(s.begin(), s.end());

		temp.clear();
		for (auto n : test){
			inData.workData.edges.push_back(n);
			temp.push_back(n);
		}
		std::cout << "\t" << i << "-Dimensional Simplices: " << temp.size() << std::endl;
		
	}
	std::cout << "\tComplex Size: " << inData.workData.edges.size() << std::endl;
	
	return inData;
}


// configPipe -> configure the function settings of this pipeline segment
bool ripsPipe::configPipe(std::map<std::string, std::string> configMap){
	
	
	auto pipe = configMap.find("dimensions");
	if(pipe != configMap.end()){
		dim = std::atoi(configMap["dimensions"].c_str());
	}
	
	
	return true;
}

