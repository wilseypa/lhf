#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <regex>
#include <fstream>
#include <cmath>
#include <map>

// Header file for readInput class - see readInput.cpp for descriptions

class readInput
{
private:
  FILE *pFile;
  char streamBuffer[1000];
  static bool parseDoubleVector(const std::string &, std::vector<double> &);

public:
  readInput();
  static std::vector<std::vector<double>> readCSV(const std::string &filename);
  static std::vector<std::vector<double>> readMAT(const std::string &filename);
  static std::map<std::vector<short>, short> readBinaryMap(const std::string &filename, int numProcesses = 1, int rank = 0);
  static std::vector<std::vector<short>> readBinaryVector(const std::string &filename);
  bool streamInit(const std::string &filename);
  bool streamRead(std::vector<double> &);
};
