#pragma once
#include "utils.hpp"

// Header file for writeOutput class - see writeOutput.cpp for descriptions

class writeOutput
{
private:
public:
  writeOutput();
  static bool writeStats(const std::string &, const std::string &);
  static bool writeRunLog(const std::string &, const std::string &);
  static bool writeCSV(const std::string &, const std::string &);
  static bool writeCSV(const std::string &, const std::string &, const std::string &);
  static bool writeCSV(std::vector<std::vector<double>> &, const std::string &);
  static bool writeCSV(std::vector<std::vector<double>> &, const std::string &, const std::string &);
  static bool writeMAT(std::vector<std::vector<double>> &, const std::string &);
  static bool writeBarcodes(std::vector<bettiBoundaryTableEntry> &, const std::string &);
  static bool writeConsole(std::vector<bettiBoundaryTableEntry> &);
  static std::string logRun(std::map<std::string, std::string> &, const std::string &, const std::string &, const std::string &);
};
