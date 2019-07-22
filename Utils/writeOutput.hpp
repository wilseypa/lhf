#pragma once

// Header file for writeOutput class - see writeOutput.cpp for descriptions

class writeOutput {
  private:
  public:
    writeOutput();
    bool writeStats(std::string, std::string);
    bool writeCSV(std::string, std::vector<std::vector<double>>);
    bool writeMAT(std::string, std::vector<std::vector<double>>);
	//bool writeConsole(pipePacket*);
};

