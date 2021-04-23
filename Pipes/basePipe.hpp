#pragma once

// Header file for basePipe class - see basePipe.cpp for descriptions
#include <map>
#include "pipePacket.hpp"
#include "utils.hpp"

class basePipe {
  private:
  public:
  	bool configured = false;
    std::string fnmod = "";
    utils ut;
    std::string pipeType = "basePipe";
    bool debug = 0;
    std::string simplicialComplex = "";
    std::string complexType = "";
    std::string outputFile;
    
    basePipe(){};
    virtual ~basePipe(){};
    static basePipe* newPipe(const std::string&, const std::string&);
    void runPipeWrapper(pipePacket&);
    virtual void outputData(pipePacket&);
    virtual void runPipe(pipePacket&);
    virtual bool configPipe(std::map<std::string, std::string>&);    
    
};

