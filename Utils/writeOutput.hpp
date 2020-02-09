#pragma once

// Header file for writeOutput class - see writeOutput.cpp for descriptions

class writeOutput {
  private:
  public:
    writeOutput();
    bool writeStats(std::string, std::string);
	bool writeCSV(std::string, std::string);
	bool writeCSV(std::string, std::string, std::string);
    bool writeCSV(std::vector<std::vector<double>>, std::string);
    bool writeCSV(std::vector<std::vector<double>>, std::string, std::string);
    bool writeMAT(std::vector<std::vector<double>>, std::string);
	bool writeBarcodes(std::vector<std::vector<double>>, std::string);
	bool writeBarcodes(std::string, std::string);
	//bool writeConsole(pipePacket*);
};

