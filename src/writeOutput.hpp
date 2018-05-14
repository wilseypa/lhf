#pragma once
#include "pipePacket.hpp"

// Header file for writeOutput class - see writeOutput.cpp for descriptions

class writeOutput {
  private:
  public:
    writeOutput();
    bool writeCSV(std::string filename, std::vector<std::vector<double>>);
    bool writeMAT(std::string filename, std::vector<std::vector<double>>);
	bool writeConsole(pipePacket*);
};

