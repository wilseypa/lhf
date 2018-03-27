#include <iostream>
#include <vector>
#include "readInput.hpp"

using namespace std;

int main(int argc, char* argv[]){
    readInput *rs = new readInput();
    vector<vector<double>> result = rs->readCSV("./testData.csv");

	cout << "_________CSV RESULTS__________" << endl;
	
	for(int i = 0; i < result.size(); i++){
		for(int j = 0; j < result[i].size(); j++){
			cout << result[i][j] << "\t";
		}
		cout << endl;
	}

	for(int i = 0; i < argc; i++){
		cout << "Arg:\t" << argv[i] << endl;
	}

    return 0;
}
