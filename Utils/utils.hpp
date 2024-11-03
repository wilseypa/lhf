#pragma once

#include <set>
#include <vector>
#include <memory>
#include <map>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <numeric>
#include "kdTree.hpp"

/**
 * @brief Simplex Node Structure
 *
 */
struct simplexNode
{
	unsigned index;
	long long hash = -1;

	std::set<unsigned> simplex = {};
	double weight = 0;
	simplexNode() {}
	simplexNode(std::set<unsigned> simp, double wt) : simplex(simp), weight(wt) {}
};

/**
 * @brief Alpha Node Structure
 *
 */
struct alphaNode
{
	unsigned index;
	long long hash = -1;

	std::set<unsigned> simplex = {};
	double weight = 0;
	double filterationvalue = 0;
	double circumRadius = 0;
	double betaRadius = 0;
	std::vector<double> betaradii;
	double volume;
	std::vector<double> hpcoff; // cofficient of simplex hyperplane
	std::vector<double> circumCenter;
	std::vector<std::vector<double>> betaCenters;
	alphaNode() {}
	alphaNode(std::set<unsigned> simp, double wt) : simplex(simp), weight(wt) {}
};

/**
 * @brief Witness Node Structure
 *
 */
struct witnessNode
{
	unsigned index;
	long long hash = -1;

	std::set<unsigned> witnessPts;
	std::vector<double> landmarkPt;

	std::set<unsigned> simplex = {};
	double weight = 0;
	witnessNode() {}
	witnessNode(std::set<unsigned> simp, double wt) : simplex(simp), weight(wt) {}
};

/**
 * @brief
 *
 */
struct bettiBoundaryTableEntry
{
	unsigned bettiDim;
	double birth;
	double death;
	std::set<unsigned> boundaryPoints;
	bool isCentroid = false;

	double getSize()
	{
		return (sizeof(double) * 2) + ((boundaryPoints.size() + 1) * sizeof(unsigned));
	}
};

class utils
{
private:
	std::string debug = "0";
	std::string outputFile = "console";

public:
	utils();
	utils(const std::string &, const std::string &);

	// Utility functions for writing to console/debug file
	void writeLog(const std::string &module, const std::string &message);
	void writeDebug(const std::string &module, const std::string &message);
	void writeError(const std::string &module, const std::string &error)
	{
		writeLog(module, error);
		return;
	};
	void writeFile(const std::string &fullMessage);

	static double computeMaxRadius(int, const std::vector<std::vector<double>> &, const std::vector<std::vector<double>> &, const std::vector<unsigned> &);
	static double computeAvgRadius(int, const std::vector<std::vector<double>> &, const std::vector<std::vector<double>> &, const std::vector<unsigned> &);
	static std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> separatePartitions(int, const std::vector<std::vector<double>> &, const std::vector<unsigned> &);
	static std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> separatePartitions(double, const std::vector<std::vector<double>> &, const std::vector<std::vector<double>> &, const std::vector<unsigned> &);
	static std::vector<std::vector<std::vector<double>>> separateBoundaryPartitions(const std::vector<std::set<unsigned>> &, const std::vector<std::vector<double>> &, const std::vector<unsigned> &);
	// void extractBoundaryPoints(std::vector<bettiBoundaryTableEntry>&);

	template <typename T>
	static std::set<unsigned> extractBoundaryPoints(const std::vector<std::shared_ptr<T>> &);

	template <typename T>
	static std::set<unsigned> extractBoundaryPoints(const std::vector<T *> &);

	static std::vector<bettiBoundaryTableEntry> mapPartitionIndexing(const std::vector<unsigned> &, std::vector<bettiBoundaryTableEntry>); // Const
	static void print2DVector(const std::vector<std::vector<unsigned>> &);
	static void print1DVector(const std::vector<unsigned> &);
	static void print1DVector(const std::set<unsigned> &);
	static void print1DVector(const std::vector<double> &);
	static void print1DSet(const std::pair<std::set<unsigned>, double> &);
	static inline double vectors_distance(const double &, const double &);
	static double vectors_distance(const std::vector<double> &, const std::vector<double> &);
	static std::set<unsigned> setXOR(const std::set<unsigned> &, const std::set<unsigned> &);
	static std::set<unsigned> setIntersect(const std::set<unsigned> &, const std::set<unsigned> &);
	static std::vector<unsigned> setIntersect(std::vector<unsigned>, std::vector<unsigned>, bool);
	static std::vector<std::set<unsigned>> getSubsets(const std::set<unsigned> &, size_t);
	static std::vector<std::set<unsigned>> getSubsets(const std::set<unsigned> &);
	static std::vector<std::vector<unsigned>> getSubsets(const std::vector<unsigned> &);
	static std::vector<unsigned> symmetricDiff(std::vector<unsigned>, std::vector<unsigned>, bool);
	static std::set<unsigned> symmetricDiff(const std::set<unsigned> &, const std::set<unsigned> &);
	static std::vector<unsigned> setUnion(std::vector<unsigned>, std::vector<unsigned>, bool);
	static std::set<unsigned> setUnion(const std::set<unsigned> &, const std::set<unsigned> &);
	static std::pair<std::vector<unsigned>, std::vector<unsigned>> intersect(std::vector<unsigned>, std::vector<unsigned>, bool);
	static bool isSubset(std::vector<unsigned>, std::vector<unsigned>);

	static bool sortBySecond(const std::pair<std::set<unsigned>, double> &, const std::pair<std::set<unsigned>, double> &);

	static double determinantOfMatrix(std::vector<std::vector<double>> mat, unsigned n); // Const
	// Alpha (delaunay)
	static double circumRadius(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> *distMatrix);
	static double circumRadius(const std::vector<short> &simplex, const std::vector<std::vector<double>> &distMatrix);
	static std::vector<double> circumCenter(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> &inputData);
	static std::vector<double> circumCenter(const std::vector<short> &simplex, const std::vector<std::vector<double>> &inputData);
	static double simplexVolume(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> *distMatrix, int dd);
	static double simplexVolume(const std::vector<std::vector<double>> &mat);
	static std::vector<std::vector<double>> inverseOfMatrix(std::vector<std::vector<double>> mat, int n); // Const
	static std::vector<std::vector<double>> matrixMultiplication(const std::vector<std::vector<double>> &matA, const std::vector<std::vector<double>> &matB);
	static std::pair<std::vector<double>, std::vector<std::vector<double>>> nullSpaceOfMatrix(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> &inputdata, std::vector<double> cc, double radius, bool lowerdimension = false); // Const

	static std::vector<std::vector<bool>> betaNeighbors(const std::vector<std::vector<double>> &, double beta, const std::string &betaMode);
	static std::vector<std::vector<double>> betaCentersCalculation(const std::vector<double> &hpcoff, double beta, double circumRadius, const std::vector<double> &circumCenter);
	static std::pair<std::vector<std::vector<double>>, std::vector<double>> calculateBetaCentersandRadius(const std::vector<unsigned> &simplex, const std::vector<std::vector<double>> &inputData, const std::vector<std::vector<double>> *distMatrix, double beta);

	static std::vector<double> serialize(const std::vector<std::vector<double>> &);
	static std::vector<std::vector<double>> deserialize(const std::vector<double> &, unsigned);

	static std::vector<double> nearestNeighbors(const std::vector<double> &, const std::vector<std::vector<double>> &);

	static Eigen::MatrixXd covariance(const std::vector<std::vector<double>> &);

	static double getAverage(const std::vector<double> &);
	static std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &);
	static std::vector<std::vector<double>> transposeMeanAdjusted(const std::vector<std::vector<double>> &);
	static std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> computePCA(const std::vector<std::vector<double>> &, int);
	static std::vector<std::vector<double>> computePCAInverse(const std::vector<std::vector<double>> &, const std::vector<std::vector<double>> &, const std::vector<std::vector<double>> &);

	template <typename T>
	static T dot(const std::vector<T> &, const std::vector<T> &);
};
