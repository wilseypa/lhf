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
  static bool writeCSV(const std::vector<std::vector<double>> &, const std::string &);
  static bool writeCSV(const std::vector<std::vector<double>> &, const std::string &, const std::string &);
  static void writeBinary(const std::map<std::vector<short>, short> &, const std::string &);
  static void writeBinary(const std::vector<std::vector<short>> &, const std::string &);
  static bool writeMAT(const std::vector<std::vector<double>> &, const std::string &);
  static bool writeBarcodes(const std::vector<bettiBoundaryTableEntry> &, const std::string &);
  static bool writeConsole(const std::vector<bettiBoundaryTableEntry> &);
  static std::string logRun(const std::map<std::string, std::string> &, const std::string &, const std::string &, const std::string &);
};
