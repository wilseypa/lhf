#include<iostream>
#include "polytopalcomplex.hpp"

using namespace std;


int main(){
	
	auto rs = readInput();
	std::vector<std::vector<double>> data;
	std::string filename = "inputfile.txt";
	//cout<<"Enter File Name::";
	//std::getline (std::cin,filename);
	data = rs.readCSV(filename);
	
	polytopalComplex polyComplex(data);
	polyComplex.persistence();
	
	for(auto x: polyComplex.getbettiTable()){
		for(auto y:x)
			cout<<y<<",";
		cout<<"\n";
	}
	return 0;
}
