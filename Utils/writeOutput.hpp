#pragma once
#include "utils.hpp"

// Header file for writeOutput class - see writeOutput.cpp for descriptions

class writeOutput {
  private:
  public:
    writeOutput();
    static bool writeStats(std::string, std::string);
    static bool writeCSV(std::string, std::string);
    static bool writeCSV(std::string, std::string, std::string);
    static bool writeCSV(std::vector<std::vector<double>>, std::string);
    static bool writeCSV(std::vector<std::vector<double>>, std::string, std::string);
    static bool writeMAT(std::vector<std::vector<double>>, std::string);
    static bool writeBarcodes(std::vector<bettiBoundaryTableEntry>, std::string);
    //bool writeConsole(pipePacket*);
};

