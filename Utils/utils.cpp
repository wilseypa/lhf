/**
 * @file utils.hpp
 *
 * @brief Contains utility functions for the LHF system (https://github.com/wilseypa/LHF).
 */

#include <string>
#include <cmath>
#include <algorithm>
#include <utility>
#include <numeric>
#include <iostream>
#include <fstream>
#include <iterator>
#include "utils.hpp"
#include <time.h>
#include <execution>
#include <Eigen/Dense>

/**
 * @brief Construct a new utils::utils object
 *
 */
utils::utils() {}

/**
 * @brief Construct a new utils::utils object
 *
 * @param _debug The debug level.
 * @param _outputFile The output file.
 */
utils::utils(const std::string &_debug, const std::string &_outputFile) : debug(_debug), outputFile(_outputFile) {}

/**
 * @brief
 *
 * @param module
 * @param message
 */
void utils::writeLog(const std::string &module, const std::string &message)
{
	if (debug == "1" || debug == "true")
	{
		std::cout << "[" << module << "]:\t" << message << std::endl;
	}
	else
	{
		writeFile("[" + module + "]:\t" + message);
	}
	return;
}

/**
 * @brief
 *
 * @param module
 * @param message
 */
void utils::writeDebug(const std::string &module, const std::string &message)
{
	if (debug == "0" || debug == "false")
	{
		return;
	}
	else
	{
		std::cout << "[DEBUG]\t[" << module << "]:\t" << message << std::endl;
		writeFile("[DEBUG]\t[" + module + "]:\t" + message);
	}

	return;
}

/**
 * @brief
 *
 * @param fullMessage
 */
void utils::writeFile(const std::string &fullMessage)
{
	std::ofstream outfile;
	outfile.open(outputFile + "_debug.txt", std::ios_base::app);
	outfile << fullMessage << "\n";

	return;
}

/**
 * @brief
 *
 * @param k
 * @param centroids
 * @param originalData
 * @param labels
 * @return double
 */
double utils::computeMaxRadius(int k, const std::vector<std::vector<double>> &centroids, const std::vector<std::vector<double>> &originalData, const std::vector<unsigned> &labels)
{
	double maxRadius = 0;
	double curRadius = 0;

	// Iterate through each point
	for (unsigned i = 0; i < originalData.size(); i++)
	{
		// Check the distance of this point to it's centroid
		curRadius = vectors_distance(originalData[i], centroids[labels[i]]);

		if (curRadius > maxRadius)
			maxRadius = curRadius;
	}

	return maxRadius;
}

/**
 * @brief
 *
 * @param k
 * @param centroids
 * @param originalData
 * @param labels
 * @return double
 */
double utils::computeAvgRadius(int k, const std::vector<std::vector<double>> &centroids, const std::vector<std::vector<double>> &originalData, const std::vector<unsigned> &labels)
{
	double totalRadius = 0;

	// Iterate through each point
	for (unsigned i = 0; i < originalData.size(); i++)
	{

		// Check the distance of this point to it's centroid
		totalRadius += vectors_distance(originalData[i], centroids[labels[i]]);
	}

	return totalRadius / originalData.size();
}

/**
 * @brief
 *
 * @param k
 * @param originalData
 * @param labels
 * @return std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>>
 */
std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> utils::separatePartitions(int k, const std::vector<std::vector<double>> &originalData, const std::vector<unsigned> &labels)
{
	std::vector<std::vector<double>> a;
	std::vector<unsigned> b;
	std::vector<std::vector<std::vector<double>>> res(k, a);
	std::vector<std::vector<unsigned>> labres(k, b);

	for (unsigned i = 0; i < labels.size(); i++)
	{
		res[labels[i]].push_back(originalData[i]);
		labres[labels[i]].push_back(i);
	}

	return std::make_pair(labres, res);
}

/**
 * @brief
 *
 * @param rad
 * @param centroids
 * @param originalData
 * @param labels
 * @return std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>>
 */
std::pair<std::vector<std::vector<unsigned>>, std::vector<std::vector<std::vector<double>>>> utils::separatePartitions(double rad, const std::vector<std::vector<double>> &centroids, const std::vector<std::vector<double>> &originalData, const std::vector<unsigned> &labels)
{
	std::vector<std::vector<double>> a;
	std::vector<unsigned> b;

	// Store the results to return
	std::vector<std::vector<std::vector<double>>> res(centroids.size(), a);
	std::vector<std::vector<unsigned>> labres(centroids.size(), b);

	// Iterate through each label
	for (unsigned i = 0; i < labels.size(); i++)
	{

		// Check for this point belonging to each centroid
		for (unsigned j = 0; j < centroids.size(); j++)
		{
			if (labels[i] == j)
			{
				// If this is a labeled constituent put to front of array
				res[j].insert(res[j].begin(), originalData[i]);
				labres[j].insert(labres[j].begin(), i);
			}
			else
			{

				double curRad = vectors_distance(originalData[i], centroids[j]);

				// If this distance is less than our radius cutoff push to back
				if (curRad < rad)
				{
					res[j].push_back(originalData[i]);
					labres[j].push_back(i);
				}
			}
		}
	}

	return std::make_pair(labres, res);
}

/**
 * @brief
 *
 * @param boundaryLists
 * @param originalData
 * @param labels
 * @return std::vector<std::vector<std::vector<double>>>
 */
std::vector<std::vector<std::vector<double>>> utils::separateBoundaryPartitions(const std::vector<std::set<unsigned>> &boundaryLists, const std::vector<std::vector<double>> &originalData, const std::vector<unsigned> &labels)
{
	std::vector<std::vector<double>> a;
	std::vector<std::vector<std::vector<double>>> res(boundaryLists.size(), a);

	for (unsigned i = 0; i < originalData.size(); i++)
	{

		for (unsigned j = 0; j < boundaryLists.size(); j++)
		{

			if (boundaryLists[j].find(labels[i]) != boundaryLists[j].end())
				res[j].push_back(originalData[i]);
		}
	}

	return res;
}

template std::set<unsigned> utils::extractBoundaryPoints<simplexNode>(const std::vector<std::shared_ptr<simplexNode>> &);
template std::set<unsigned> utils::extractBoundaryPoints<alphaNode>(const std::vector<std::shared_ptr<alphaNode>> &);
template std::set<unsigned> utils::extractBoundaryPoints<witnessNode>(const std::vector<std::shared_ptr<witnessNode>> &);

/**
 * @brief
 *
 * @tparam T
 * @param boundary
 * @return std::set<unsigned>
 */
template <typename T>
std::set<unsigned> utils::extractBoundaryPoints(const std::vector<std::shared_ptr<T>> &boundary)
{
	std::set<unsigned> boundaryPoints;
	for (auto &simplex : boundary)
		boundaryPoints.insert(simplex->simplex.begin(), simplex->simplex.end());
	return boundaryPoints;
}

template std::set<unsigned> utils::extractBoundaryPoints<simplexNode>(const std::vector<simplexNode *> &);
template std::set<unsigned> utils::extractBoundaryPoints<alphaNode>(const std::vector<alphaNode *> &);
template std::set<unsigned> utils::extractBoundaryPoints<witnessNode>(const std::vector<witnessNode *> &);

/**
 * @brief
 *
 * @tparam T
 * @param boundary
 * @return std::set<unsigned>
 */
template <typename T>
std::set<unsigned> utils::extractBoundaryPoints(const std::vector<T *> &boundary)
{
	std::set<unsigned> boundaryPoints;
	for (auto &simplex : boundary)
		boundaryPoints.insert(simplex->simplex.begin(), simplex->simplex.end());
	return boundaryPoints;
}

/**
 * @brief
 *
 * @param v The vector of doubles.
 * @return double The average of the vector of doubles.
 */
double utils ::getAverage(const std::vector<double> &v)
{
	if (v.empty())
	{
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

/**
 * @brief
 *
 * @param input
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> utils ::transposeMeanAdjusted(const std::vector<std::vector<double>> &input)
{
	std::vector<double> inputtranspose(input.size());
	std::vector<std::vector<double>> inputstd(input[0].size(), std::vector<double>(input.size()));

	for (size_t i = 0; i < input[0].size(); i++)
	{
		for (size_t j = 0; j < input.size(); j++)
			inputtranspose[j] = input[j][i];
		double average = getAverage(inputtranspose);
		for (size_t j = 0; j < inputtranspose.size(); j++)
			inputstd[i][j] = inputtranspose[j] - average;
	}
	return inputstd;
}

/**
 * @brief
 *
 * @param input
 * @return Eigen::MatrixXd
 */
Eigen::MatrixXd utils ::covariance(const std::vector<std::vector<double>> &input)
{
	Eigen::MatrixXd result(input.size(), input.size());
	int i = 0, j = 0;
	for (auto &x : input)
	{
		for (auto &y : input)
		{
			double value = 0;
			double averagex = getAverage(x);
			double averagey = getAverage(y);
			for (size_t i = 0; i < x.size(); i++)
			{
				value += (x[i] - averagex) * (y[i] - averagey);
			}
			result(i, j) = value / (x.size() - 1);
			j++;
		}
		i++;
		j = 0;
	}
	return result;
}

/**
 * @brief
 *
 * @param mat
 * @param n
 * @return double
 */
double utils ::determinantOfMatrix(std::vector<std::vector<double>> mat, unsigned n)
{
	double det = 1;
	unsigned index;
	for (unsigned i = 0; i < n; i++)
	{
		index = i;
		while (mat[index][i] == 0 && index < n)
			index++;
		if (index == n)
			continue;
		if (index != i)
		{
			for (unsigned j = 0; j < n; j++)
			{
				double temp12 = mat[index][j];
				mat[index][j] = mat[i][j];
				mat[i][j] = temp12;
			}
			det = det * pow(-1, index - i);
		}
		double rectemp = mat[i][i];
		for (unsigned j = i; j < n; j++)
			mat[i][j] /= rectemp;
		for (unsigned j = i + 1; j < n; j++)
		{
			if (mat[j][i] != 0)
			{
				double rectemp2 = mat[j][i];
				for (unsigned t = i; t < n; t++)
					mat[i][t] *= rectemp2;
				for (unsigned t = i; t < n; t++)
					mat[j][t] -= mat[i][t];
				for (unsigned k = 0; k < n; k++)
					mat[i][k] /= rectemp2;
			}
		}
		for (unsigned k = 0; k < n; k++)
			mat[i][k] *= rectemp;
	}
	for (unsigned i = 0; i < n; i++)
		det = det * mat[i][i];

	return (det);
}

/**
 * @brief
 *
 * @param input
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> utils ::transpose(const std::vector<std::vector<double>> &input)
{
	std::vector<std::vector<double>> inputtranspose(input[0].size(), std::vector<double>(input.size()));
	for (size_t i = 0; i < input[0].size(); i++)
	{
		for (size_t j = 0; j < input.size(); j++)
			inputtranspose[i][j] = input[j][i];
	}
	return inputtranspose;
}

/**
 * @brief
 *
 * @param matA
 * @param matB
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> utils ::matrixMultiplication(const std::vector<std::vector<double>> &matA, const std::vector<std::vector<double>> &matB)
{
	size_t n1 = matA.size();
	size_t m1 = matA[0].size();
	size_t n2 = matB.size();
	size_t m2 = matB[0].size();
	std::vector<std::vector<double>> mat(n1, std::vector<double>(m2, 0));

	if (m1 != n2)
		return mat;

	for (size_t i = 0; i < n1; i++)
	{
		for (size_t j = 0; j < m2; j++)
		{
			for (size_t x = 0; x < m1; x++)
			{
				mat[i][j] += matA[i][x] * matB[x][j];
			}
		}
	}
	return mat;
}

/**
 * @brief
 *
 * @param mat
 * @param n
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> utils ::inverseOfMatrix(std::vector<std::vector<double>> mat, int n)
{
	int index;
	std::vector<std::vector<double>> matinv(n, std::vector<double>(n, 0));

	for (int i = 0; i < n; i++)
		matinv[i][i] = 1;

	for (int i = 0; i < n; i++)
	{
		index = i;
		while (mat[index][i] == 0 && index < n)
			index++;
		if (index == n)
			continue;
		if (index != i)
		{
			for (int j = 0; j < n; j++)
			{
				double temp12 = mat[index][j];
				mat[index][j] = mat[i][j];
				mat[i][j] = temp12;
				double temp121 = matinv[index][j];
				matinv[index][j] = matinv[i][j];
				matinv[i][j] = temp121;
			}
		}
		double rectemp = mat[i][i];
		if (mat[i][i] != 1)
		{
			for (int j = 0; j < n; j++)
			{
				mat[i][j] /= rectemp;
				matinv[i][j] /= rectemp;
			}
		}
		for (int j = 0; j < n; j++)
		{
			if (mat[j][i] != 0 && j != i)
			{
				double rectemp2 = mat[j][i];
				for (int t = 0; t < n; t++)
				{
					mat[i][t] *= rectemp2;
					matinv[i][t] *= rectemp2;
				}
				for (int t = 0; t < n; t++)
				{
					mat[j][t] -= mat[i][t];
					matinv[j][t] -= matinv[i][t];
				}
				for (int k = 0; k < n; k++)
				{
					mat[i][k] /= rectemp2;
					matinv[i][k] /= rectemp2;
				}
			}
		}
	}
	return matinv;
}

/**
 * @brief
 *
 * @param hpcoff
 * @param beta
 * @param circumRadius
 * @param circumCenter
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> utils ::betaCentersCalculation(const std::vector<double> &hpcoff, double beta, double circumRadius, const std::vector<double> &circumCenter)
{
	double distance = sqrt(pow((beta * circumRadius), 2) - pow(circumRadius, 2));
	double d1, d2; // Parallel Plane coefficient
	double sqrtofsquaredsum = 0, squaredsum = 0;
	double dotproduct = 0;
	int i = 0;
	for (auto x : hpcoff)
	{
		squaredsum += x * x;
		dotproduct += x * circumCenter[i];
		i++;
	}

	sqrtofsquaredsum = sqrt(squaredsum);

	d1 = -dotproduct + distance * sqrtofsquaredsum;
	d2 = -dotproduct - distance * sqrtofsquaredsum;

	double t1, t2;

	t1 = (-dotproduct - d1) / squaredsum;

	t2 = (-dotproduct - d2) / squaredsum;

	std::vector<std::vector<double>> centers;
	std::vector<double> center1;
	std::vector<double> center2;
	i = 0;
	for (auto x : hpcoff)
	{
		center1.push_back(x * t1 + circumCenter[i]);
		center2.push_back(x * t2 + circumCenter[i]);
		i++;
	}
	centers.push_back(center1);
	centers.push_back(center2);
	return centers;
}

/**
 * @brief
 *
 * @param input
 * @param targetDim
 * @return std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
 */
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> utils::computePCA(const std::vector<std::vector<double>> &input, int targetDim)
{
	std::vector<std::vector<double>> transp = transposeMeanAdjusted(input);
	Eigen::MatrixXd covar = covariance(transp);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
	es.compute(covar);
	std::vector<double> eigenval;
	std::vector<std::vector<double>> eigenvec;
	for (int k = input[0].size() - 1; k >= 0; k--)
	{
		std::vector<double> evec;
		for (size_t p = 0; p < input[0].size(); p++)
			evec.push_back(es.eigenvectors().col(k)[p]);
		eigenvec.push_back(evec);
	}

	for (int p = input[0].size() - 1; p >= 0; p--)
		eigenval.push_back(es.eigenvalues()[p]);
	std::vector<std::vector<double>> eigenVecSelect;

	for (int k = 0; k < targetDim; k++)
		eigenVecSelect.push_back(eigenvec[k]);

	return std::make_pair(transpose(matrixMultiplication(eigenVecSelect, transp)), eigenVecSelect);
}

/**
 * @brief
 *
 * @param input
 * @param FinalOutput
 * @param eigenvectors
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> utils::computePCAInverse(const std::vector<std::vector<double>> &input, const std::vector<std::vector<double>> &FinalOutput, const std::vector<std::vector<double>> &eigenvectors)
{

	auto result = transpose(matrixMultiplication(transpose(eigenvectors), transpose(FinalOutput)));

	std::vector<std::vector<double>> transposeinput = transpose(input);
	int i = 0;
	for (auto x : transposeinput)
	{
		double average = getAverage(x);
		for (size_t y = 0; y < result.size(); y++)
		{
			result[y][i] += average;
		}
		i++;
	}
	return result;
}

/**
 * @brief
 *
 * @param simplex
 * @param inputData
 * @param cc
 * @param radius
 * @param lowerdimension
 * @return std::pair<std::vector<double>, std::vector<std::vector<double>>>
 */
std::pair<std::vector<double>, std::vector<std::vector<double>>> utils ::nullSpaceOfMatrix(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> &inputData, std::vector<double> cc, double radius, bool lowerdimension)
{
	int index;
	srand(time(NULL));
	int n = simplex.size();
	std::vector<double> matns(n, 1);
	std::vector<std::vector<double>> mat;

	for (auto x : simplex)
		mat.push_back(inputData[x]);
	mat.push_back(cc);
	std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> outputPCA;
	if (lowerdimension == true)
	{
		outputPCA = computePCA(mat, simplex.size());
		mat = outputPCA.first;
		cc = mat[mat.size() - 1];
		mat.pop_back();
	}

	for (int i = 0; i < n; i++)
	{
		index = i;
		while (mat[index][i] == 0 && index < n)
			index++;
		if (index == n)
			continue;
		if (index != i)
		{
			for (int j = 0; j < n; j++)
			{
				double temp12 = mat[index][j];
				mat[index][j] = mat[i][j];
				mat[i][j] = temp12;
			}
			double temp121 = matns[index];
			matns[index] = matns[i];
			matns[i] = temp121;
		}
		double rectemp = mat[i][i];
		if (mat[i][i] != 1)
		{
			for (int j = 0; j < n; j++)
			{
				mat[i][j] /= rectemp;
			}
			matns[i] /= rectemp;
		}

		for (int j = 0; j < n; j++)
		{
			if (mat[j][i] != 0 && j != i)
			{
				double rectemp2 = mat[j][i];
				for (int t = 0; t < n; t++)
				{
					mat[i][t] *= rectemp2;
				}
				matns[i] *= rectemp2;
				for (int t = 0; t < n; t++)
				{
					mat[j][t] -= mat[i][t];
				}
				matns[j] -= matns[i];
				for (int k = 0; k < n; k++)
				{
					mat[i][k] /= rectemp2;
				}
				matns[i] /= rectemp2;
			}
		}
	}

	return std::make_pair(matns, outputPCA.second);
}

/**
 * @brief
 *
 * @param simplex
 * @param inputData
 * @return std::vector<double>
 */
std::vector<double> utils::circumCenter(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> &inputData)
{
	const int n = inputData[0].size();
	const int m = simplex.size();

	Eigen::MatrixXd matA(m, m);
	Eigen::VectorXd matC(m, 1);
	std::vector<double> circumCenter(n, 0.0);

	std::set<unsigned> simplexcopy = simplex;

	auto it = simplexcopy.end();
	it--;
	const unsigned Sn = *(it);
	simplexcopy.erase(Sn);

	int ii = 0;
	for (auto i : simplexcopy)
	{
		int temp = ii;
		const Eigen::VectorXd d1 = Eigen::Map<const Eigen::VectorXd>(inputData[i].data(), n) - Eigen::Map<const Eigen::VectorXd>(inputData[Sn].data(), n);
		for (auto j : simplexcopy)
		{
			if (i <= j)
			{
				const Eigen::VectorXd d2 = Eigen::Map<const Eigen::VectorXd>(inputData[j].data(), n) - Eigen::Map<const Eigen::VectorXd>(inputData[Sn].data(), n);
				const double dotProduct = d1.dot(d2);
				matA(ii, temp) = dotProduct;
				matA(temp, ii) = dotProduct;

				if (i == j)
					matC(ii) = dotProduct / 2.0;

				temp++;
			}
		}
		matA(ii, m - 1) = 0;
		ii++;
	}

	matA.row(ii).setConstant(1);
	matC(ii) = 1;

	// Solve the linear system
	Eigen::VectorXd rawCircumCenter = matA.inverse() * matC;

	for (int i = 0; i < n; i++)
	{
		auto index = simplex.begin();
		for (int j = 0; j < m; j++, ++index)
		{
			circumCenter[i] += rawCircumCenter(j) * inputData[*index][i];
		}
	}

	return circumCenter;
}

/**
 * @brief
 *
 * @param simplex
 * @param distMatrix
 * @return double
 */
double utils::circumRadius(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> *distMatrix)
{
	unsigned n = simplex.size();
	Eigen::MatrixXd matA(n, n);
	Eigen::MatrixXd matACap(n + 1, n + 1);
	unsigned ii = 0;
	for (auto i : simplex)
	{
		matACap.row(ii + 1).col(0).setConstant(1); // Set column 0 of matACap to 1
		unsigned temp = 0;
		for (auto j : simplex)
		{
			double distSquared = pow(((*distMatrix)[i][j] != 0) ? (*distMatrix)[i][j] : (*distMatrix)[j][i], 2);
			matA(ii, temp) = distSquared;
			matACap(ii + 1, temp + 1) = distSquared;
			temp++;
		}
		ii++;
	}
	matACap.row(0).setConstant(1);
	matACap(0, 0) = 0;
	double result = (-matA.determinant() / (2 * matACap.determinant()));
	if (result < 0)
		result = std::numeric_limits<double>::max();
	return result;
}

/**
 * @brief
 *
 * @param simplex
 * @param inputData
 * @return std::vector<double>
 */
std::vector<double> utils::circumCenter(const std::vector<short> &simplex, const std::vector<std::vector<double>> &inputData)
{
	const int n = inputData[0].size();
	const int m = simplex.size();
	Eigen::MatrixXd matA(m, m);
	Eigen::VectorXd matC(m, 1);
	std::vector<double> circumCenter(n, 0.0);
	const Eigen::VectorXd &Sn = Eigen::Map<const Eigen::VectorXd>(inputData[simplex.back()].data(), n);
	for (int i = 0; i < m - 1; ++i)
	{
		const Eigen::VectorXd &d1 = Eigen::Map<const Eigen::VectorXd>(inputData[simplex[i]].data(), n) - Sn;
		for (int j = i; j < m - 1; ++j)
		{
			const double dotProduct = d1.dot(Eigen::Map<const Eigen::VectorXd>(inputData[simplex[j]].data(), n) - Sn);
			matA(i, j) = dotProduct;
			matA(j, i) = dotProduct;
		}
	}
	matA.col(m - 1).setZero();
	matA.row(m - 1).setOnes();
	matC = matA.diagonal() / 2;
	matC(m - 1) = 1;
	// Solve the linear
	Eigen::VectorXd rawCircumCenter = matA.inverse() * matC;
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			circumCenter[j] += rawCircumCenter(i) * inputData[simplex[i]][j];
	return circumCenter;
}

/**
 * @brief
 *
 * @param simplex
 * @param distMatrix
 * @return double
 */
double utils::circumRadius(const std::vector<short> &simplex, const std::vector<std::vector<double>> &distMatrix)
{
	const unsigned n = simplex.size();
	Eigen::MatrixXd matA(n, n);
	Eigen::MatrixXd matACap(n + 1, n + 1);
	for (unsigned i = 0; i < n; ++i)
		for (unsigned j = 0; j < n; ++j)
			matA(i, j) = pow(simplex[i] <= simplex[j] ? distMatrix[simplex[i]][simplex[j]] : distMatrix[simplex[j]][simplex[i]], 2);
	matACap.block(1, 1, n, n) = matA;
	matACap.col(0).setConstant(1);
	matACap.row(0).setConstant(1);
	matACap(0, 0) = 0;
	double result = -matA.determinant() / (2 * matACap.determinant());
	if (result < 0)
	{
		std::numeric_limits<double>::max();
	}
	return result;
}

/**
 * @brief
 *
 * @param simplex
 * @param distMatrix
 * @param dd
 * @return double
 */
double utils ::simplexVolume(const std::set<unsigned> &simplex, const std::vector<std::vector<double>> *distMatrix, int dd)
{
	std::vector<std::vector<double>> matACap(simplex.size() + 1);
	int ii = 0;
	for (auto i : simplex)
	{
		matACap[ii + 1].push_back(1);
		for (auto j : simplex)
		{
			if ((*distMatrix)[i][j] != 0)
			{
				matACap[ii + 1].push_back(pow(((*distMatrix)[i][j]), 2));
			}
			else
			{
				matACap[ii + 1].push_back(pow(((*distMatrix)[j][i]), 2));
			}
		}
		ii++;
	}
	matACap[0].push_back(0);
	for (size_t i = 0; i < simplex.size(); i++)
		matACap[0].push_back(1);
	if (dd % 2 == 0)
		return ((-1) * determinantOfMatrix(matACap, simplex.size() + 1)) / (pow(2, dd) * (pow(tgamma(dd + 1), 2)));
	else
		return (determinantOfMatrix(matACap, simplex.size() + 1)) / (pow(2, dd) * (pow(tgamma(dd + 1), 2)));
}

/**
 * @brief
 *
 * @param spoints
 * @return double
 */
double utils ::simplexVolume(const std::vector<std::vector<double>> &spoints)
{
	std::vector<std::vector<double>> matACap(spoints.size() + 1);
	int ii = 0;
	for (auto i : spoints)
	{
		matACap[ii + 1].push_back(1);
		for (auto j : spoints)
		{
			matACap[ii + 1].push_back(pow(vectors_distance(i, j), 2));
		}
		ii++;
	}
	matACap[0].push_back(0);
	for (size_t i = 0; i < spoints.size(); i++)
		matACap[0].push_back(1);
	if (spoints.size() % 2 == 0)
		return ((-1) * determinantOfMatrix(matACap, spoints.size() + 1)) / (pow(2, spoints[0].size()) * (pow(tgamma(spoints[0].size() + 1), 2)));
	else
		return (determinantOfMatrix(matACap, spoints.size() + 1)) / (pow(2, spoints[0].size()) * (pow(tgamma(spoints[0].size() + 1), 2)));
}

/**
 * @brief
 *
 * @param inData
 * @param beta
 * @param betaMode
 * @return std::vector<std::vector<bool>>
 */
std::vector<std::vector<bool>> utils ::betaNeighbors(const std::vector<std::vector<double>> &inData, double beta, const std::string &betaMode)
{
	std::vector<std::vector<bool>> incidenceMatrix(inData.size(), std::vector<bool>(inData.size(), 0));
	kdTree tree(inData, inData.size()); // KDTree for efficient nearest neighbor search
	double betafactor;
	std::vector<std::vector<double>> rotationMatrix(inData[0].size(), std::vector<double>(inData[0].size(), 0));

	//[
	// cos(90) -sin(90)  ................... 0
	// sin(90) cos(90)  .....................0
	// 0		0	0 ......0.......0
	//......................1......	0	0
	// 0		0	0	1	0
	// 0		0	0	0    	1
	//]

	for (unsigned i = 0; i < inData[0].size(); i++)
	{
		for (unsigned j = 0; j < inData[0].size(); j++)
		{
			if (i == 1 && j == 0)
				rotationMatrix[i][j] = 1;
			else if (i == 0 && j == 1)
				rotationMatrix[i][j] = -1;
			else if (i == 0 && j == 0)
				rotationMatrix[i][j] = 0;
			else if (i == 1 && j == 1)
				rotationMatrix[i][j] = 0;
			else if (i == j)
				rotationMatrix[i][j] = 1;
			else
				rotationMatrix[i][j] = 0;
		}
	}

	if (beta < 0)
	{
		std::cout << "Invalid Beta Value :: Value set to 1" << std::endl;
		beta = 1; // Changing to Default
	}
	for (unsigned i = 0; i < inData.size(); i++)
	{
		for (unsigned j = i + 1; j < inData.size(); j++)
		{
			/*	std::cout<<" First Point \n";
				for(auto x:inData[i])
					std::cout<<x<<" ";

				std::cout<<" Second Point \n";
				for(auto x:inData[j])
					std::cout<<x<<" ";
			*/
			double radius;
			std::vector<double> temp;
			std::vector<double> center1;
			std::vector<double> center2;
			std::vector<double> temp1;
			std::vector<double> tempfinal;
			std::vector<std::vector<double>> colvector;
			if (beta <= 1)
			{ // Common Definition for lune and circle based beta-Skeleton
				radius = vectors_distance(inData[i], inData[j]) / (2 * beta);
				std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp), [](double e1, double e2)
							   { return ((e1 + e2) / 2); });
				betafactor = sqrt((1 - pow(beta, 2))) / (2 * beta);
				std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp1), [](double e1, double e2)
							   { return ((e2 - e1)); });
				for (auto x : temp1)
				{
					std::vector<double> singlevector;
					singlevector.push_back(x);
					colvector.push_back(singlevector);
				}
				std::vector<std::vector<double>> matmul = matrixMultiplication(rotationMatrix, colvector);
				for (auto x : matmul)
					for (auto y : x)
						tempfinal.push_back(y * betafactor);
				std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center1), [](double e1, double e2)
							   { return ((e1 + e2)); });
				std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center2), [](double e1, double e2)
							   { return ((e1 - e2)); });

				/*	     	   std::cout<<" First Center and radius:: "<<radius<<"\n";
						   for(auto x:center1)
							std::cout<<x<<" ";
							   std::cout<<" Second Center \n";
						   for(auto x:center2)
							std::cout<<x<<" ";
				*/
				std::vector<size_t> neighbors1 = tree.neighborhoodIndices(center1, radius); // All neighbors in radius-ball
				std::vector<unsigned> toremove;
				toremove.push_back(i);
				toremove.push_back(j);
				std::sort(toremove.begin(), toremove.end());
				std::sort(neighbors1.begin(), neighbors1.end());
				std::vector<int> neg1;
				std::set_difference(neighbors1.begin(), neighbors1.end(), toremove.begin(), toremove.end(), std::back_inserter(neg1));

				std::vector<size_t> neighbors2 = tree.neighborhoodIndices(center2, radius); // All neighbors in radius-ball
				std::sort(neighbors2.begin(), neighbors2.end());
				std::vector<int> neg2;
				std::set_difference(neighbors2.begin(), neighbors2.end(), toremove.begin(), toremove.end(), std::back_inserter(neg2));

				std::vector<size_t> neighbors;

				std::vector<size_t> v(std::min(neg1.size(), neg2.size()));
				std::vector<size_t>::iterator it;
				std::sort(neg1.begin(), neg1.end());
				std::sort(neg2.begin(), neg2.end());
				it = std::set_intersection(neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
				v.resize(it - v.begin());
				neighbors = v;

				if (neighbors.size() == 0)
					incidenceMatrix[i][j] = true;
				else
					incidenceMatrix[i][j] = false;
			}
			else if (betaMode == "lune")
			{ // Lune Based Beta Skeleton for beta > 1
				radius = (vectors_distance(inData[i], inData[j]) * beta) / 2;
				double bf = beta / 2;
				double bf1 = 1 - beta / 2;

				std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(center1), [bf, bf1](double e1, double e2)
							   { return ((e1 * bf + e2 * bf1)); });
				std::transform(inData[j].begin(), inData[j].end(), inData[i].begin(), std::back_inserter(center2), [bf, bf1](double e1, double e2)
							   { return ((e1 * bf + e2 * bf1)); });

				std::vector<size_t> neighbors1 = tree.neighborhoodIndices(center1, radius); // All neighbors in radius-ball
				std::vector<unsigned> toremove;
				toremove.push_back(i);
				toremove.push_back(j);
				std::sort(toremove.begin(), toremove.end());
				std::sort(neighbors1.begin(), neighbors1.end());
				std::vector<int> neg1;
				std::set_difference(neighbors1.begin(), neighbors1.end(), toremove.begin(), toremove.end(), std::back_inserter(neg1));

				std::vector<size_t> neighbors2 = tree.neighborhoodIndices(center2, radius); // All neighbors in radius-ball
				std::sort(neighbors2.begin(), neighbors2.end());
				std::vector<int> neg2;
				std::set_difference(neighbors2.begin(), neighbors2.end(), toremove.begin(), toremove.end(), std::back_inserter(neg2));

				std::vector<size_t> neighbors;

				std::vector<size_t> v(std::min(neg1.size(), neg2.size()));
				std::vector<size_t>::iterator it;
				std::sort(neg1.begin(), neg1.end());
				std::sort(neg2.begin(), neg2.end());
				it = std::set_intersection(neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
				v.resize(it - v.begin());
				neighbors = v;

				if (neighbors.size() == 0)
					incidenceMatrix[i][j] = true;
				else
					incidenceMatrix[i][j] = false;
			}
			else if (betaMode == "circle")
			{ // Circle Based Beta Skeleton for beta >1
				radius = (vectors_distance(inData[i], inData[j]) * beta) / (2);
				std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp), [](double e1, double e2)
							   { return ((e1 + e2) / 2); });
				betafactor = sqrt((pow(beta, 2) - 1)) / 2;
				std::transform(inData[i].begin(), inData[i].end(), inData[j].begin(), std::back_inserter(temp1), [](double e1, double e2)
							   { return ((e2 - e1)); });
				for (auto x : temp1)
				{
					std::vector<double> singlevector;
					singlevector.push_back(x);
					colvector.push_back(singlevector);
				}
				std::vector<std::vector<double>> matmul = matrixMultiplication(rotationMatrix, colvector);
				for (auto x : matmul)
					for (auto y : x)
						tempfinal.push_back(y * betafactor);
				std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center1), [](double e1, double e2)
							   { return ((e1 + e2)); });
				std::transform(temp.begin(), temp.end(), tempfinal.begin(), std::back_inserter(center2), [](double e1, double e2)
							   { return ((e1 - e2)); });

				std::vector<size_t> neighbors1 = tree.neighborhoodIndices(center1, radius); // All neighbors in radius-ball
				std::vector<unsigned> toremove;
				toremove.push_back(i);
				toremove.push_back(j);
				std::sort(toremove.begin(), toremove.end());
				std::sort(neighbors1.begin(), neighbors1.end());
				std::vector<int> neg1;
				std::set_difference(neighbors1.begin(), neighbors1.end(), toremove.begin(), toremove.end(), std::back_inserter(neg1));

				std::vector<size_t> neighbors2 = tree.neighborhoodIndices(center2, radius); // All neighbors in radius-ball
				std::sort(neighbors2.begin(), neighbors2.end());
				std::vector<int> neg2;
				std::set_difference(neighbors2.begin(), neighbors2.end(), toremove.begin(), toremove.end(), std::back_inserter(neg2));

				std::vector<size_t> neighbors;

				std::vector<size_t> v(neighbors1.size() + neighbors2.size());
				std::vector<size_t>::iterator it;
				std::sort(neg1.begin(), neg1.end());
				std::sort(neg2.begin(), neg2.end());
				it = std::set_union(neg1.begin(), neg1.end(), neg2.begin(), neg2.end(), v.begin());
				v.resize(it - v.begin());
				neighbors = v;
				if (neighbors.size() == 0)
					incidenceMatrix[i][j] = true;
				else
					incidenceMatrix[i][j] = false;
			}
		}
	}
	return incidenceMatrix;
}

// void utils::extractBoundaryPoints(std::vector<bettiBoundaryTableEntry>& bettiTable){
// 	for(auto& bet: bettiTable){
// 		std::set<unsigned> bound;
// 		for(auto simplex : bet.boundary) bound.insert(simplex->simplex.begin(), simplex->simplex.end());
// 		bet.boundaryPoints = bound;
// 	}
// }

/**
 * @brief
 *
 * @param partitionedLabels
 * @param bettiTable
 * @return std::vector<bettiBoundaryTableEntry>
 */
std::vector<bettiBoundaryTableEntry> utils::mapPartitionIndexing(const std::vector<unsigned> &partitionedLabels, std::vector<bettiBoundaryTableEntry> bettiTable)
{
	for (auto &bet : bettiTable)
	{
		std::set<unsigned> convBound;

		for (auto ind : bet.boundaryPoints)
		{
			convBound.insert(partitionedLabels[ind]);
		}

		bet.boundaryPoints = convBound;
	}
	return bettiTable;
}

/**
 * @brief
 *
 * @param a
 */
void utils::print2DVector(const std::vector<std::vector<unsigned>> &a)
{
	for (unsigned i = 0; i < a.size(); i++)
	{
		for (unsigned j = 0; j < a[i].size(); j++)
		{
			std::cout << a[i][j] << '\t';
		}
		std::cout << std::endl;
	}
	return;
}

/**
 * @brief
 *
 * @param a
 */
void utils::print1DSet(const std::pair<std::set<unsigned>, double> &a)
{
	std::cout << "Test\t";

	for (auto iter = a.first.begin(); iter != a.first.end(); iter++)
	{
		std::cout << *iter << ",";
	}
	std::cout << "\t";
}

/**
 * @brief
 *
 * @param a
 */
void utils::print1DVector(const std::vector<double> &a)
{
	for (unsigned i = 0; i < a.size(); i++)
	{
		std::cout << a[i] << ",";
	}
	std::cout << "\n";
	return;
}

/**
 * @brief
 *
 * @param a
 */
void utils::print1DVector(const std::vector<unsigned> &a)
{
	for (unsigned i = 0; i < a.size(); i++)
	{
		std::cout << a[i] << ",";
	}
	std::cout << "\n";
	return;
}

/**
 * @brief
 *
 * @param a
 */
void utils::print1DVector(const std::set<unsigned> &a)
{
	for (auto &z : a)
	{
		std::cout << z << ",";
	}
	std::cout << "\n";
	return;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return double
 */
inline double utils::vectors_distance(const double &a, const double &b)
{
	return pow((a - b), 2);
}

/**
 * @brief
 *
 * @param setA
 * @param setB
 * @return std::set<unsigned>
 */
std::set<unsigned> utils::setXOR(const std::set<unsigned> &setA, const std::set<unsigned> &setB)
{
	std::set<unsigned> ret;

	std::set_symmetric_difference(setA.begin(), setA.end(), setB.begin(), setB.end(), std::inserter(ret, ret.begin()));

	return ret;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return double
 */
double utils::vectors_distance(const std::vector<double> &a, const std::vector<double> &b)
{
	if (b.size() == 0)
		return 0;

#ifndef NO_PARALLEL_ALGORITHMS
	return sqrt(std::transform_reduce(std::execution::par, a.cbegin(), a.cend(), b.cbegin(), 0.0, std::plus<>(),
									  [](double e1, double e2)
									  { return (e1 - e2) * (e1 - e2); }));
#else
	return sqrt(std::transform_reduce(a.cbegin(), a.cend(), b.cbegin(), 0.0, std::plus<>(),
									  [](double e1, double e2)
									  { return (e1 - e2) * (e1 - e2); }));
#endif
}

/**
 * @brief
 *
 * @param v1
 * @param v2
 * @param isSorted
 * @return std::vector<unsigned>
 */
std::vector<unsigned> utils::setIntersect(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted)
{
	std::vector<unsigned> ret;

	if (v1 == v2)
		return v1;

	if (!isSorted)
	{
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}

	std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(ret));

	return ret;
}

/**
 * @brief
 *
 * @param v1
 * @param v2
 * @param isSorted
 * @return std::set<unsigned>
 */
std::set<unsigned> utils::setIntersect(const std::set<unsigned> &v1, const std::set<unsigned> &v2)
{
	std::set<unsigned> ret;

	if (v1 == v2)
		return v1;

	std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(ret, ret.begin()));

	return ret;
}

/**
 * @brief
 *
 * @param v1
 * @param v2
 * @param isSorted
 * @return std::pair<std::vector<unsigned>, std::vector<unsigned>>
 */
// Find the intersect of two vectors
std::pair<std::vector<unsigned>, std::vector<unsigned>> utils::intersect(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted)
{
	std::pair<std::vector<unsigned>, std::vector<unsigned>> ret;

	if (v1 == v2)
		return ret;

	if (!isSorted)
	{
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}

	std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(ret.first));

	if (ret.first.size() == 2)
	{
		std::set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(ret.second));
		return ret;
	}
	else
		return std::pair<std::vector<unsigned>, std::vector<unsigned>>();
}

/**
 * @brief
 *
 * @param v1
 * @param v2
 * @param isSorted
 * @return std::vector<unsigned>
 */
// Find the symmetric difference of two vectors
std::vector<unsigned> utils::symmetricDiff(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted)
{
	std::vector<unsigned> ret;
	std::vector<unsigned> retTemp;

	if (v1 == v2)
		return ret;

	if (!isSorted)
	{
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}
	std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(retTemp));

	return retTemp;
}

/**
 * @brief
 *
 * @param v1
 * @param v2
 * @param isSorted
 * @return std::set<unsigned>
 */
// Find the symmetric difference of two vectors
std::set<unsigned> utils::symmetricDiff(const std::set<unsigned> &v1, const std::set<unsigned> &v2)
{
	std::set<unsigned> ret;
	if (v1 == v2)
		return ret;

	std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(ret, ret.begin()));
	return ret;
}

/**
 * @brief
 *
 * @param set
 * @param dim
 * @return std::vector<std::set<unsigned>>
 */
// Iteratively build subsets (faces) of the simplex set
std::vector<std::set<unsigned>> utils::getSubsets(const std::set<unsigned> &set, size_t dim)
{
	std::vector<std::set<unsigned>> subset;
	subset.push_back(std::set<unsigned>());

	// For each set in the
	for (auto i = set.begin(); i != set.end(); i++)
	{
		std::vector<std::set<unsigned>> subsetTemp = subset;

		for (unsigned j = 0; j < subsetTemp.size(); j++)
		{
			subsetTemp[j].insert(*i);
		}

		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++)
		{
			subset.push_back(*j);
		}
	}

	std::vector<std::set<unsigned>> retSubset;

	for (std::set<unsigned> z : subset)
	{
		if (z.size() == dim)
			retSubset.push_back(z);
	}
	return retSubset;
}

/**
 * @brief
 *
 * @param v1
 * @param v2
 * @return std::set<unsigned>
 */
// Find the union of two vectors
std::set<unsigned> utils::setUnion(const std::set<unsigned> &v1, const std::set<unsigned> &v2)
{
	std::set<unsigned> retTemp;

	if (v1 == v2)
		return v1;

	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), std::inserter(retTemp, retTemp.begin()));

	return retTemp;
}

/**
 * @brief
 *
 * @param v1
 * @param v2
 * @param isSorted
 * @return std::vector<unsigned>
 */
// Find the union of two vectors
std::vector<unsigned> utils::setUnion(std::vector<unsigned> v1, std::vector<unsigned> v2, bool isSorted)
{
	std::vector<unsigned> ret;
	std::vector<unsigned> retTemp;

	if (v1 == v2)
		return ret;

	if (!isSorted)
	{
		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());
	}

	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(retTemp));

	/*for(auto iter = v1.begin(); iter!= v1.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = v2.begin(); iter!= v2.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << "\t";
	for(auto iter = retTemp.begin(); iter!= retTemp.end(); iter++){
		std::cout << *iter << ",";
	}
	std::cout << std::endl;*/
	return retTemp;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return true
 * @return false
 */
bool utils::sortBySecond(const std::pair<std::set<unsigned>, double> &a, const std::pair<std::set<unsigned>, double> &b)
{
	if (a.second != b.second)
		return a.second < b.second;
	// If the simplices have the same weight, sort by reverse lexicographic order for fastPersistence
	auto itA = a.first.rbegin(), itB = b.first.rbegin();
	while (itA != a.first.rend())
	{
		if (*itA != *itB)
			return *itA > *itB;
		++itA;
		++itB;
	}
	return false;
}

/**
 * @brief
 *
 * @param set
 * @return std::vector<std::set<unsigned>>
 */
// Iteratively build subsets (faces) of the simplex set
std::vector<std::set<unsigned>> utils::getSubsets(const std::set<unsigned> &set)
{
	std::vector<std::set<unsigned>> subset;
	std::set<unsigned> empty;
	subset.push_back(empty);

	// For each set in the
	for (auto i = set.begin(); i != set.end(); i++)
	{
		std::vector<std::set<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++)
		{
			subsetTemp[j].insert(entry);
		}

		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++)
		{
			subset.push_back(*j);
		}
	}

	std::vector<std::set<unsigned>> retSubset;

	for (std::set<unsigned> z : subset)
	{
		if (z.size() == set.size() - 1)
			retSubset.push_back(z);
	}
	return retSubset;
}

/**
 * @brief
 *
 * @param set
 * @return std::vector<std::vector<unsigned>>
 */
// Iteratively build subsets (faces) of the simplex set
std::vector<std::vector<unsigned>> utils::getSubsets(const std::vector<unsigned> &set)
{
	std::vector<std::vector<unsigned>> subset;
	std::vector<unsigned> empty;
	subset.push_back(empty);

	// For each set in the
	for (auto i = set.begin(); i != set.end(); i++)
	{
		std::vector<std::vector<unsigned>> subsetTemp = subset;
		unsigned entry = *i;

		for (unsigned j = 0; j < subsetTemp.size(); j++)
		{
			subsetTemp[j].push_back(entry);
		}

		for (auto j = subsetTemp.begin(); j != subsetTemp.end(); j++)
		{
			subset.push_back(*j);
		}
	}

	std::vector<std::vector<unsigned>> retSubset;

	for (std::vector<unsigned> z : subset)
	{
		if (z.size() == set.size() - 1)
			retSubset.push_back(z);
	}
	return retSubset;
}

/**
 * @brief
 *
 * @param point
 * @param pointcloud
 * @return std::vector<double>
 */
std::vector<double> utils::nearestNeighbors(const std::vector<double> &point, const std::vector<std::vector<double>> &pointcloud)
{
	// based on random projection, x is current point being examined, n is number of centroids/facilities
	std::vector<double> retVal(pointcloud.size());

	// Get sq distances for each point
	for (auto &currentPoint : pointcloud)
	{
		retVal.push_back(vectors_distance(point, currentPoint));
	}

	return retVal;
}

/**
 * @brief
 *
 * @param serialData
 * @param dim
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>> utils::deserialize(const std::vector<double> &serialData, unsigned dim)
{

	// First check if the vector size matches the dimension
	if (serialData.size() % dim != 0)
	{
		std::cout << "Error occurred when deserializing data: invalid size" << std::endl;
		return {};
	}

	// Deduce the number of vectors
	size_t n = serialData.size() / dim;

	std::vector<std::vector<double>> ret(n, std::vector<double>(dim));

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			ret[i][j] = serialData[(i * dim) + j];
		}
	}
	return ret;
}

/**
 * @brief
 *
 * @param origData
 * @return std::vector<double>
 */
std::vector<double> utils::serialize(const std::vector<std::vector<double>> &origData)
{

	// Make sure we have data to serialize
	if (origData.size() == 0)
	{
		std::cout << "Error occurred when serializing data: empty vector argument" << std::endl;
		return {};
	}

	size_t n = origData.size();
	size_t d = origData[0].size();

	// Preallocate our vector to prevent resizing
	std::vector<double> ret(n * d);

	// Copy element by element
	for (size_t i = 0; i < n; i++)
	{
		for (size_t k = 0; k < d; k++)
		{
			ret[(i * d) + k] = origData[i][k];
		}
	}

	return ret;
}

/**
 * @brief
 *
 * @param dsimplex
 * @param inputData
 * @param distMatrix
 * @param beta
 * @return std::pair<std::vector<std::vector<double>>, std::vector<double>>
 */
std::pair<std::vector<std::vector<double>>, std::vector<double>> utils::calculateBetaCentersandRadius(const std::vector<unsigned> &dsimplex, const std::vector<std::vector<double>> &inputData, const std::vector<std::vector<double>> *distMatrix, double beta)
{
	std::vector<std::vector<double>> betacenters;
	std::vector<double> betaradii;
	bool intersectionCircle = false;
	if (beta < 0)
		exit(0);
	else if (beta == 0)
		return std::make_pair(betacenters, betaradii);
	else if (beta < 1)
		intersectionCircle = true;

	if (beta < 1)
		beta = 1 / beta;
	std::set<unsigned> simplex(dsimplex.begin(), dsimplex.end());
	std::vector<double> circumCenter;
	if (simplex.size() > 2)
		circumCenter = utils::circumCenter(simplex, inputData);
	else if (simplex.size() == 2)
	{
		auto first = simplex.begin();
		std::vector<double> R;
		std::vector<double> A = inputData[*first];
		std::advance(first, 1);
		std::vector<double> B = inputData[*first];
		std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(R), [](double e1, double e2)
					   { return ((e1 + e2) / 2); });
		circumCenter = R;
	}

	double circumRadius;
	if (simplex.size() > 2)
		circumRadius = utils::circumRadius(simplex, distMatrix);
	else
		circumRadius = pow((*distMatrix)[dsimplex[0]][dsimplex[1]] / 2, 2);

	std::vector<size_t> neighbors;
	std::vector<std::vector<size_t>> neighborsCircleIntersection;
	for (auto x : simplex)
	{
		double expr1, expr2, expr3;
		std::vector<unsigned> face1;
		face1 = dsimplex;
		face1.erase(std::remove(face1.begin(), face1.end(), x), face1.end());
		std::set<unsigned> face(face1.begin(), face1.end());
		std::vector<double> faceCC;
		if (face.size() > 2)
			faceCC = utils::circumCenter(face, inputData);
		double faceRadius;
		if (face.size() > 2)
			faceRadius = utils::circumRadius(face, distMatrix);
		else
			faceRadius = pow((*distMatrix)[face1[0]][face1[1]] / 2, 2);
		auto result = utils::nullSpaceOfMatrix(face, inputData, faceCC, sqrt(faceRadius));
		std::vector<double> hpcoff = result.first;
		std::vector<double> betaCenter;
		double betaRadius;
		std::vector<std::vector<double>> betaCenters;
		bool sameside = false;
		if (intersectionCircle && beta > 2)
		{
			if (beta < 3)
			{
				double ratio = sqrt(circumRadius) / sqrt(faceRadius);
				betaRadius = sqrt(faceRadius) + (beta - 2) * (sqrt(circumRadius) - sqrt(faceRadius));
				betaCenters = utils::betaCentersCalculation(hpcoff, 1 + (beta - 2) * (ratio - 1), sqrt(faceRadius), faceCC);
			}
			else
			{
				betaCenters = utils::betaCentersCalculation(hpcoff, beta - 1, sqrt(faceRadius), faceCC);
				betaRadius = sqrt(faceRadius) * (beta - 1);
			}
			expr1 = 0;
			expr2 = 0;
			expr3 = 0;
			for (unsigned i = 0; i < hpcoff.size(); i++)
			{
				expr1 += hpcoff[i] * circumCenter[i];
				expr2 += hpcoff[i] * betaCenters[0][i];
				expr3 += hpcoff[i] * inputData[x][i];
			}
			expr1--;
			expr2--;
			expr3--;
			if (((expr1 > 0 && expr3 > 0) && expr2 > 0) || ((expr1 < 0 && expr3 < 0) && expr2 < 0))
			{
				sameside = true;
				betaCenter = betaCenters[1];
			}
			else if ((expr1 > 0 && expr3 > 0) || (expr1 < 0 && expr3 < 0))
			{
				sameside = true;
				if (expr2 > 0 && expr1 > 0)
					betaCenter = betaCenters[1];
				else
					betaCenter = betaCenters[0];
			}
			else
			{
				sameside = false;
				if (expr2 > 0 && expr1 > 0)
					betaCenter = betaCenters[0];
				else
					betaCenter = betaCenters[1];
			}
		}
		else
		{
			expr1 = 0;
			expr3 = 0;
			for (unsigned i = 0; i < hpcoff.size(); i++)
			{
				expr1 += hpcoff[i] * circumCenter[i];
				expr3 += hpcoff[i] * inputData[x][i];
			}
			expr1--;
			expr3--;
			if ((expr1 > 0 && expr3 > 0) || (expr1 < 0 && expr3 < 0))
				sameside = true;
			else
				sameside = false;
			for (unsigned y = 0; y < inputData[0].size(); y++)
				if (sameside)
				{
					if (intersectionCircle)
						betaCenter.push_back((2 - beta) * circumCenter[y] + (beta - 1) * faceCC[y]);
					else
						betaCenter.push_back(beta * circumCenter[y] - (beta - 1) * faceCC[y]);
				}
				else
				{
					if (intersectionCircle)
					{
						betaCenter.push_back(beta * circumCenter[y] - (beta - 1) * faceCC[y]);
					}
					else
						betaCenter.push_back((2 - beta) * circumCenter[y] + (beta - 1) * faceCC[y]);
				}
			// Anurag -> Author:: Beta Radius might remain unintialized here
		}
		if (!intersectionCircle || beta <= 2)
			betaRadius = utils::vectors_distance(betaCenter, inputData[face1[0]]);

		betacenters.push_back(betaCenter);
		betaradii.push_back(betaRadius);
	}
	return std::make_pair(betacenters, betaradii);
}

/**
 * @brief
 *
 * @param A
 * @param B
 * @return true
 * @return false
 */
bool utils::isSubset(std::vector<unsigned> A, std::vector<unsigned> B)
{
	std::sort(A.begin(), A.end());
	std::sort(B.begin(), B.end());
	return std::includes(A.begin(), A.end(), B.begin(), B.end());
}

template <typename T>
inline T utils::dot(const std::vector<T> &a, const std::vector<T> &b) { return std::inner_product(a.begin(), a.end(), b.begin(), static_cast<T>(0)); }

template double utils::dot<double>(const std::vector<double> &a, const std::vector<double> &);