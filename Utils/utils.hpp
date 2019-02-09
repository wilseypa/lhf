#pragma once

// Header file for utils class - see utils.cpp for descriptions

class utils {
  private:
  public:
	utils();
	void print2DVector(const std::vector<std::vector<unsigned>>&);
	void print1DVector(const std::vector<unsigned>&);
	double vectors_distance(const std::vector<double>&, const std::vector<double>&);
	void print1DSet(const auto&);
	std::vector<unsigned> symmetricDiff(std::vector<unsigned>, std::vector<unsigned>, bool);
	std::vector<unsigned> setUnion(std::vector<unsigned>, std::vector<unsigned>, bool);
	std::pair<std::vector<unsigned>, std::vector<unsigned>> intersect(std::vector<unsigned>, std::vector<unsigned>, bool);
};

