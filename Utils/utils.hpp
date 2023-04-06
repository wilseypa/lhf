#pragma once

#include <set>
#include <vector>
#include <memory>
#include <map>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Eigenvalues>
#include <cmath>
#include <numeric>
#include "../Preprocessing/kdTree.hpp"


// Simplex Node Structure
struct simplexNode{
	unsigned index;
	long long hash = -1;
	
	std::set<unsigned> simplex = {};
	double weight = 0;
	simplexNode(){}
	simplexNode(std::set<unsigned> simp, double wt) : simplex(simp), weight(wt) {}
};

// Alpha Node Structure
struct alphaNode{
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
	alphaNode(){}
	alphaNode(std::set<unsigned> simp, double wt) : simplex(simp), weight(wt) {}
};

// Witness Node Structure
struct witnessNode{
	unsigned index;
	long long hash = -1;
	
	std::set<unsigned> witnessPts;
	std::vector<double> landmarkPt;
	
	std::set<unsigned> simplex = {};
	double weight = 0;
	witnessNode(){}
	witnessNode(std::set<unsigned> simp, double wt) : simplex(simp), weight(wt) {}
};

struct bettiBoundaryTableEntry{
	unsigned bettiDim;
	double birth;
	double death;
	std::set<unsigned> boundaryPoints;
	bool isCentroid = false;

	double getSize(){
		return (sizeof(double)*2) + ((boundaryPoints.size() + 1) * sizeof(unsigned));
	}
};


class unionFind{
	private:
		std::vector<int> rank, parent;
	public:
		unionFind(int n);
		int find(int i);
		bool join(int x, int y);
};

class utils {
  private:
	std::string debug = "0";
	std::string outputFile = "console";

  public:
	utils();
	utils(std::string, std::string);

	static double computeMaxRadius(int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<unsigned>&);
	static double computeAvgRadius(int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<unsigned>&);
	static std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> separatePartitions(int, std::vector<std::vector<double>>, std::vector<unsigned>);
	static std::vector<std::vector<std::vector<double>>> separateBoundaryPartitions(std::vector<std::set<unsigned>>, std::vector<std::vector<double>>, std::vector<unsigned>);
	static std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> separatePartitions(double, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<unsigned>);
	// void extractBoundaryPoints(std::vector<bettiBoundaryTableEntry>&);
	
	template <typename T>
	std::set<unsigned> extractBoundaryPoints(std::vector<std::shared_ptr<T>>);
	
	template <typename T>
	std::set<unsigned> extractBoundaryPoints(std::vector<T*>);
	
	static std::vector<bettiBoundaryTableEntry> mapPartitionIndexing(std::vector<unsigned>, std::vector<bettiBoundaryTableEntry>);
	static void print2DVector(const std::vector<std::vector<unsigned>>&);
	static void print1DVector(const std::vector<unsigned>&);
	static void print1DVector(const std::set<unsigned>&);
	static void print1DVector(const std::vector<double>&);
	static double vectors_distance(const double&, const double&);
	static double vectors_distance(const std::vector<double>&, const std::vector<double>&);
	static void print1DSet(const std::pair<std::set<unsigned>, double>&);
	static std::set<unsigned> setXOR(std::set<unsigned>&, std::set<unsigned>&);
	static std::set<unsigned> setIntersect(std::set<unsigned>, std::set<unsigned>, bool isSorted);
	static std::vector<unsigned> setIntersect(std::vector<unsigned>, std::vector<unsigned>, bool);
	static std::vector<std::set<unsigned>> getSubsets(std::set<unsigned>, int);
	static std::vector<unsigned> symmetricDiff(std::vector<unsigned>, std::vector<unsigned>, bool);
	static std::set<unsigned> symmetricDiff(std::set<unsigned>, std::set<unsigned>, bool);
	static std::vector<unsigned> setUnion(std::vector<unsigned>, std::vector<unsigned>, bool);
	static std::set<unsigned> setUnion(std::set<unsigned>, std::set<unsigned>);
	static std::pair<std::vector<unsigned>, std::vector<unsigned>> intersect(std::vector<unsigned>, std::vector<unsigned>, bool);
	static bool isSubset(std::vector<unsigned>,std::vector<unsigned>);
	//Utility functions for writing to console/debug file
	void writeLog(std::string module, std::string message);
	void writeDebug(std::string module, std::string message);
	void writeError(std::string module, std::string error){writeLog(module,error);return;};
	void writeFile(std::string fullMessage);

	static bool sortBySecond(const std::pair<std::set<unsigned>, double> &, const std::pair<std::set<unsigned>, double> &);
	static std::vector<std::set<unsigned>> getSubsets(std::set<unsigned> set);
	static std::vector<std::vector<unsigned>> getSubsets(std::vector<unsigned> set);

	static double determinantOfMatrix(std::vector<std::vector<double>> mat, int n);
	
	//Alpha (delaunay)
	static double circumRadius(std::set<unsigned> simplex,std::vector<std::vector<double>>* distMatrix);
	static std::vector<double> circumCenter(std::set<unsigned> simplex,std::vector<std::vector<double>> inputData);
	
	
	static std::vector<std::vector<double>> inverseOfMatrix(std::vector<std::vector<double>> mat, int n);
	static std::vector<std::vector<double>> matrixMultiplication(std::vector<std::vector<double>> matA, std::vector<std::vector<double>> matB);
    static std::pair<std::vector<double>,std::vector<std::vector<double>>> nullSpaceOfMatrix(std::set<unsigned> simplex, std::vector<std::vector<double>> inputdata,std::vector<double> &cc, double radius,bool lowerdimension=false);
    static std::vector<std::vector<double>> betaCentersCalculation(std::vector<double> hpcoff, double beta, double circumRadius,std::vector<double> circumCenter);
	static double simplexVolume(std::set<unsigned> simplex,std::vector<std::vector<double>>* distMatrix,int dd);
	static double simplexVolume(std::vector<std::vector<double>>mat);
    
    static std::vector<std::vector<bool>> betaNeighbors(std::vector<std::vector<double>> &,double beta,std::string betaMode);
	static std::pair<std::vector<std::vector<double>>,std::vector<double>> calculateBetaCentersandRadius(std::vector<unsigned> simplex ,std::vector<std::vector<double>> &inputData,std::vector<std::vector<double>>* distMatrix, double beta);

	static std::vector<double> serialize(std::vector<std::vector<double>>& );
	static std::vector<std::vector<double>> deserialize(std::vector<double> , unsigned);

	static std::vector<double> nearestNeighbors(std::vector<double>&, std::vector<std::vector<double>>&);
	
	static Eigen::MatrixXd  covariance(std::vector<std::vector<double>> );
	static std::vector<std::vector<double>> transposeMeanAdjusted(std::vector<std::vector<double>> );
	static double getAverage(std::vector<double> &);
	static std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> );
    static std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>> computePCA(std::vector<std::vector<double>> , int);
    static std::vector<std::vector<double>> computePCAInverse(std::vector<std::vector<double>>,std::vector<std::vector<double>> , std::vector<std::vector<double>> );



};
