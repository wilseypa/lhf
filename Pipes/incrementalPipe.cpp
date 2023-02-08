#include "incrementalPipe.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include <fstream>
#include <omp.h>

template <typename nodeType>
int incrementalPipe<nodeType>::expand_d_minus_1_simplex(std::vector<unsigned> simplex, unsigned omission)
{
	double rad_vector[this->inputData.size()];
	std::set<unsigned> simp(simplex.begin(), simplex.end());
	simp.erase(omission);
	std::cout << "New Simplex ";
	for (auto i : simp)
		std::cout << i << " ";
	std::cout << std::endl;
	double largest_radius = 0, smallest_radius = 10000000000;
	int partitioning_point = -1, triangulation_point = -1;
	for (auto new_point : this->search_space)
	{
		if (simp.find(new_point) != simp.end())
			continue;
		simp.insert(new_point);
		rad_vector[new_point] = sqrt(utils::circumRadius(simp, this->distMatrix));
		// std::cout << new_point << " " << rad_vector[new_point] << std::endl;
		simp.erase(new_point);
		if (largest_radius < rad_vector[new_point])
		{
			largest_radius = rad_vector[new_point];
			partitioning_point = new_point;
		}
	}
	simp.insert(partitioning_point);
	std::vector<double> largest_circle_center = utils::circumCenter(simp, this->inputData);
	simp.erase(partitioning_point);

	if (largest_radius > utils::vectors_distance(largest_circle_center, this->inputData[omission])) // prior point lied inside the circle new_triangulation will be outside the largest_circle
	{
		for (auto new_point : this->search_space)
		{
			if (simp.find(new_point) == simp.end() && new_point != omission && smallest_radius > rad_vector[new_point] && largest_radius < utils::vectors_distance(largest_circle_center, this->inputData[new_point]))
			{
				smallest_radius = rad_vector[new_point];
				triangulation_point = new_point;
			}
		}
	}
	else
	{
		for (auto new_point : this->search_space)
		{
			if (simp.find(new_point) == simp.end() && new_point != omission && smallest_radius > rad_vector[new_point] && largest_radius > utils::vectors_distance(largest_circle_center, this->inputData[new_point]))
			{
				smallest_radius = rad_vector[new_point];
				triangulation_point = new_point;
			}
		}
	}
	if (triangulation_point != -1)
	{
		simp.insert(triangulation_point);
		std::vector<double> smallest_circle_center = utils::circumCenter(simp, this->inputData);
		for (auto new_point : this->search_space)
		{
			if (smallest_radius < rad_vector[new_point] && smallest_radius > utils::vectors_distance(smallest_circle_center, this->inputData[new_point]))
			{
				smallest_radius = rad_vector[new_point];
				triangulation_point = (int)new_point;
			}
		}
	}
	std::cout << triangulation_point << std::endl;
	return -1;
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
	this->inputData = inData.inputData;
	unsigned dim = inData.inputData[0].size();
	unsigned data_set_size = inData.inputData.size();
	this->search_space = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	this->distMatrix = (((alphaComplex<alphaNode> *)inData.complex)->distMatrix);
	std::vector<std::vector<unsigned>> dsimplexes;
	/*
3 4 0
3 5 4
7 3 8
9 7 8
8 3 0
6 8 0
2 9 8
6 1 8
2 8 1
*/
	std::vector<std::vector<unsigned>> inner_dsimplexes_shell = {{3, 4, 0}, {3, 5, 4}, {7, 3, 8}, {9, 7, 8}, {8, 3, 0}, {6, 8, 0}, {2, 9, 8}, {6, 1, 8}, {2, 8, 1}};
	dsimplexes.push_back(inner_dsimplexes_shell[0]);
	std::vector<std::vector<unsigned>> outer_dsimplexes_shell;
	std::vector<unsigned> simplex;
	int new_point;
	unsigned omission;
	while (inner_dsimplexes_shell.size() != 0)
	{
		outer_dsimplexes_shell.clear();
		for (int i = 0; i < inner_dsimplexes_shell.size(); i++)
		{
			for (unsigned j = 0; j <= dim; j++)
			{
				simplex = inner_dsimplexes_shell[i];
				omission = simplex[j];
				new_point = expand_d_minus_1_simplex(simplex, omission);
				if (new_point == -1)
					continue;
				simplex.erase(simplex.begin() + j);
				simplex.push_back(new_point);
				outer_dsimplexes_shell.push_back(simplex);
			}
		}
		inner_dsimplexes_shell = outer_dsimplexes_shell;
		dsimplexes.insert(dsimplexes.end(), inner_dsimplexes_shell.begin(), inner_dsimplexes_shell.end());
	}
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
	this->ut.writeDebug("delaunayPipe", "Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

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