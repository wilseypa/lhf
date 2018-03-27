#pragma once

// Header file for readInput class - see readInput.cpp for descriptions

class readInput {
  private:
  public:
    readInput();
    std::vector<std::vector<double>> readCSV(std::string filename);
    std::vector<std::vector<double>> readMAT(std::string filename);
};

