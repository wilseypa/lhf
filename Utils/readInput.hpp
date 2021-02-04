#pragma once


#include <string>
#include <iostream>
#include <vector>
#include <regex>
#include <fstream>
#include <cmath>

// Header file for readInput class - see readInput.cpp for descriptions

class readInput {
  private:
  	FILE *pFile;
  	char streamBuffer[1000];
  	static bool parseDoubleVector(std::string, std::vector<double> &);
  
  public:
    readInput();
    static std::vector<std::vector<double>> readCSV(std::string filename);
    static std::vector<std::vector<double>> readMAT(std::string filename);
    bool streamInit(std::string filename);
    bool streamRead(std::vector<double>&);
};

