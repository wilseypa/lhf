#pragma once
#include "utils.hpp"

// Header file for writeOutput class - see writeOutput.cpp for descriptions

class writeOutput {
  private:
  public:
    writeOutput();
    static bool writeStats(std::string, std::string);
    static bool writeRunLog(std::string, std::string);
    static bool writeCSV(std::string, std::string);
    static bool writeCSV(std::string, std::string, std::string);
    static bool writeCSV(std::vector<std::vector<double>>, std::string);
    static bool writeCSV(std::vector<std::vector<double>>, std::string, std::string);
    static bool writeMAT(std::vector<std::vector<double>>, std::string);
    static bool writeBarcodes(std::vector<bettiBoundaryTableEntry>, std::string);
    static bool writeConsole(std::vector<bettiBoundaryTableEntry>);
    static std::string logRun(std::map<std::string, std::string>, std::string, std::string, std::string);
};

