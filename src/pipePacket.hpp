#pragma once

// Header file for pipePacket class - see pipePacket.cpp for descriptions

class pipePacket {
  private:
  public:
	pipePacket();
  
	struct pipeData{
		std::vector<std::vector<double>> originalData;
		std::vector<std::vector<double>> workingData;
		std::vector<std::vector<int>> edges;
		std::vector<int> weights;
	} workData;
	
};

