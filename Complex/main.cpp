#include<iostream>
#include "polytopalcomplex.hpp"

using namespace std;


int main(){
	
	auto rs = readInput();
	std::vector<std::vector<double>> data;
	std::string filename = "AMSfull.csv";
	cout<<"Enter File Name::";
	std::getline (std::cin,filename);
	data = rs.readCSV(filename);
	double approxfactor = 0;
	//std::cin>>approxfactor;
	polytopalComplex polyComplex(data,approxfactor);
	polyComplex.persistenceMatrix();
	std::cout<<polyComplex.getbettiTable().size()<<" => betti Table Size\n";
	for(auto x: polyComplex.getbettiTable()){
		for(auto y:x)
			cout<<y<<",";
		cout<<"\n";
	}
	return 0;
}
