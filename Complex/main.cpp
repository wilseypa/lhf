#include<iostream>
#include "polytopalcomplex.hpp"

using namespace std;


int main(){
	
	auto rs = readInput();
	std::vector<std::vector<double>> data;
	std::string filename = "inputfile.txt";
	//std::getline (std::cin,filename);
	data = rs.readCSV(filename);
	
	polytopalComplex polyComplex(data);
	return 0;
}
