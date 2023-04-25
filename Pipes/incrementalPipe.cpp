
#include "utils.hpp"
#include <fstream>
#include <unistd.h>
#include "readInput.hpp"
#include <limits>
#include <omp.h>
#include <incrementalPipe.hpp>

template <typename T>
void printVector(const std::vector<T> &vec)
{
	std::cout << "[ ";
	for (const auto &elem : vec)
	{
		std::cout << elem << ",";
	}
	std::cout << "]\n";
}

template <typename T>
void printVector(const std::set<T> &vec)
{
	std::cout << "[ ";
	for (const auto &elem : vec)
	{
		std::cout << elem << ",";
	}
	std::cout << "]\n";
}

std::vector<unsigned> first_simplex(std::vector<std::vector<double>> &inputData, std::vector<std::vector<double>> &distMatrix)
{
	unsigned dim = inputData[0].size();
	unsigned n_pts = inputData.size();
	std::vector<unsigned> simplex;
	if (n_pts < dim + 2)
		return simplex;
	simplex.push_back(0);
	while (simplex.size() != dim + 2)
	{
		double curr_min_dist = std::numeric_limits<double>::max();
		unsigned min_dist_idx = 0;
		for (unsigned i = 1; i < n_pts; i++)
		{
			if (std::find(simplex.begin(), simplex.end(), i) != simplex.end())
				continue;
			double dist = distMatrix[0][i];
			if (dist < curr_min_dist)
			{
				curr_min_dist = dist;
				min_dist_idx = i;
			}
		}
		simplex.push_back(min_dist_idx);
	}
	sort(simplex.begin(), simplex.end());
	unsigned int pow_set_size = pow(2, simplex.size());
	std::vector<unsigned> gensimp;
	for (int counter = 1; counter < pow_set_size; counter++)
	{
		if (__builtin_popcount(counter) != dim)
			continue;
		for (int j = 0; j < simplex.size(); j++)
		{
			if (counter & (1 << j))
			{
				gensimp.push_back(simplex[j]);
			}
		}
		double radius = 0;
		std::vector<double> center;
		std::set<unsigned> simplex_set(gensimp.begin(), gensimp.end());
		for (unsigned i = 0; i < inputData.size(); i++)
		{
			if (simplex_set.find(i) != simplex_set.end())
				continue;
			simplex_set.insert(i);
			center = utils::circumCenter(simplex_set, inputData);
			radius = utils::vectors_distance(center, inputData[i]);
			for (unsigned point = 0; point < inputData.size(); point++)
			{
				if (simplex_set.find(point) != simplex_set.end())
					continue;
				if (utils::vectors_distance(center, inputData[point]) < radius)
					break;
				else if (point == inputData.size() - 1)
				{
					printVector(simplex_set);
					std::vector<unsigned> temp(simplex_set.begin(), simplex_set.end());
					return temp;
				}
			}
			simplex_set.erase(i);
		}
		gensimp.clear();
	}
	simplex.clear();
	return simplex;
}

template <typename T>
std::set<T> operator+(std::set<T> a, T b)
{
	a.insert(b);
	return a;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
	std::vector<T> temp;
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp), [](double e1, double e2)
				   { return (e1 - e2); });
	return temp;
}

double dot(const std::vector<double> &a, const std::vector<double> &b)
{
	std::vector<double> temp;
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp), [](double e1, double e2)
				   { return (e1 * e2); });
	return std::accumulate(temp.begin(), temp.end(), 0.0);
}

std::vector<std::vector<double>> distMat(std::vector<std::vector<double>> &inData)
{
	std::vector<std::vector<double>> distMatrix;
	if (distMatrix.size() > 0)
		distMatrix.clear();
	distMatrix.resize(inData.size(), std::vector<double>(inData.size(), 0));
	for (unsigned i = 0; i < inData.size(); i++)
	{
		for (unsigned j = i + 1; j < inData.size(); j++)
		{
			distMatrix[i][j] = utils::vectors_distance(inData[i], inData[j]);
		}
	}
	return distMatrix;
}

int expand_d_minus_1_simplex(std::vector<unsigned> &simp_vector, unsigned &omission, std::vector<std::vector<double>> &inputData, std::vector<unsigned> &search_space, std::vector<std::vector<double>> &distMatrix)
{
	std::set<unsigned> simp(simp_vector.begin(), simp_vector.end());
	auto p1 = utils::circumCenter(simp, inputData);
	simp.insert(omission);
	auto normal = utils::circumCenter(simp, inputData) - p1;
	simp.erase(omission);
	auto direction = (dot(normal, inputData[omission] - p1) > 0);
	double smallest_radius = std::numeric_limits<double>::max(), largest_radius = 0, ring_radius = utils::vectors_distance(p1, inputData[*(simp.begin())]);
	double curr_radius = 0;
	int triangulation_point = -1;
	bool flag = true;
	for (auto &new_point : search_space)
	{
		if (simp.find(new_point) == simp.end() && new_point != omission && !(direction ^ dot(normal, inputData[new_point] - p1) < 0))
		{
			if (ring_radius > utils::vectors_distance(p1, inputData[new_point]))
			{
				simp.insert(new_point);
				curr_radius = sqrt(utils::circumRadius(simp, &distMatrix));
				simp.erase(new_point);
				if (largest_radius < curr_radius)
				{
					largest_radius = curr_radius;
					triangulation_point = new_point;
					flag = false;
				}
			}
			else if (flag)
			{
				simp.insert(new_point);
				curr_radius = sqrt(utils::circumRadius(simp, &distMatrix));
				simp.erase(new_point);
				if (smallest_radius > curr_radius)
				{
					smallest_radius = curr_radius;
					triangulation_point = new_point;
				}
			}
		}
	}
	return triangulation_point;
}

void reduce(std::vector<std::vector<unsigned>> &outer_dsimplexes, std::vector<std::pair<std::vector<unsigned>, unsigned>> &inner_d_1_shell)
{
	std::map<std::vector<unsigned>, unsigned> outer_d_1_shell;
	for (auto &new_simplex : outer_dsimplexes)
	{
		for (auto &i : new_simplex)
		{
			std::vector<unsigned> key = new_simplex;
			key.erase(std::find(key.begin(), key.end(), i));
			if (!outer_d_1_shell.emplace(key, i).second)
				outer_d_1_shell.erase(key); // Create new shell and remove collided faces max only 2 can occur.
		}
	}
	for (auto &simp : inner_d_1_shell)
		outer_d_1_shell.erase(simp.first); // Remove faces from previous iteration
	inner_d_1_shell.clear();
	for (auto &simp : outer_d_1_shell)
		inner_d_1_shell.push_back(simp); // Remove faces from previous iteration
	outer_dsimplexes.clear();
	return;
}

// basePipe constructor
template <typename nodeType>
incrementalPipe<nodeType>::incrementalPipe()
{
  this->pipeType = "incrementalPipe";
  return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void incrementalPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	std::vector<std::vector<double>> inputData = inData.inputData;
	unsigned dim = inputData[0].size();
	unsigned data_set_size = inputData.size();
	std::vector<unsigned> search_space;
	for (unsigned i = 0; i < inputData.size(); i++)
		search_space.push_back(i);
	std::vector<std::vector<double>> distMatrix = distMat(inputData);
	std::vector<std::vector<unsigned>> dsimplexes = {first_simplex(inputData, distMatrix)}; // Get initial delaunay triangle
	std::vector<std::pair<std::vector<unsigned>, unsigned>> inner_d_1_shell;
	std::vector<std::vector<unsigned>> outer_dsimplexes;
	for (auto &new_simplex : dsimplexes)
	{
		for (auto &i : new_simplex)
		{
			std::vector<unsigned> key = new_simplex;
			key.erase(std::find(key.begin(), key.end(), i));
			inner_d_1_shell.push_back(std::make_pair(key, i));
		}
	}
	int new_point;
	while (inner_d_1_shell.size() != 0)
	{
		std::cout << dsimplexes.size() << std::endl;
#pragma omp parallel for
		for (int i = 0; i < inner_d_1_shell.size(); i++)
		{
			auto iter = inner_d_1_shell[i];
			new_point = expand_d_minus_1_simplex(iter.first, iter.second, inputData, search_space, distMatrix);
			if (new_point == -1)
				continue;
			iter.first.push_back(new_point);
			std::sort(iter.first.begin(), iter.first.end());
#pragma omp critical
			{
				if (std::find(dsimplexes.begin(), dsimplexes.end(), iter.first) == dsimplexes.end())
				{
					outer_dsimplexes.push_back(iter.first);
					dsimplexes.push_back(iter.first);
				}
			}
		}
		reduce(outer_dsimplexes, inner_d_1_shell);
	}
	std::cout << dsimplexes.size();
	return;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool incrementalPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
{
	std::string strDebug;

	auto pipe = configMap.find("debug");
	if (pipe != configMap.end())
	{
		this->debug = std::atoi(configMap["debug"].c_str());
		strDebug = configMap["debug"];
	}
	pipe = configMap.find("outputFile");
	if (pipe != configMap.end())
		this->outputFile = configMap["outputFile"].c_str();

	this->ut = utils(strDebug, this->outputFile);

	this->configured = true;
	this->ut.writeDebug("incrementalPipe", "Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void incrementalPipe<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	std::ofstream file;
	file.open("output/" + this->pipeType + "_output.csv");
	// code to print the data

	file.close();
	return;
}

template class incrementalPipe<simplexNode>;
template class incrementalPipe<alphaNode>;
template class incrementalPipe<witnessNode>;