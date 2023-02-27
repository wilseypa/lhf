#include "incrementalPipe.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include <fstream>
#include <omp.h>
#include <unistd.h>

const double tolerance = 1;

template <typename nodeType>
double incrementalPipe<nodeType>::circumRadius(std::set<unsigned> simp, unsigned new_point)
{
	simp.insert(new_point);
	return sqrt(utils::circumRadius(simp, this->distMatrix));
}

template <typename nodeType>
int incrementalPipe<nodeType>::expand_d_minus_1_simplex(std::set<unsigned> simp, unsigned omission)
{
	double rad_vector[this->inputData.size()];
	double largest_radius = -100, smallest_radius = 1e+308, ring_radius;
	int partitioning_point = -1, triangulation_point = -1;
	for (auto new_point : this->search_space)
	{
		if (simp.find(new_point) != simp.end())
			continue;
		rad_vector[new_point] = circumRadius(simp, new_point);
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
			if (simp.find(new_point) == simp.end() && new_point != omission && smallest_radius > rad_vector[new_point] && largest_radius < tolerance * utils::vectors_distance(largest_circle_center, this->inputData[new_point]))
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
			if (simp.find(new_point) == simp.end() && new_point != omission && smallest_radius > rad_vector[new_point] && largest_radius * tolerance > utils::vectors_distance(largest_circle_center, this->inputData[new_point]))
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
		ring_radius = smallest_radius;
		for (auto new_point : this->search_space)
		{
			if (simp.find(new_point) == simp.end() && new_point != omission && smallest_radius < rad_vector[new_point] && ring_radius * 0.99999999 > utils::vectors_distance(smallest_circle_center, this->inputData[new_point]))
			{
				smallest_radius = rad_vector[new_point];
				triangulation_point = (int)new_point;
			}
		}
	}
	return triangulation_point;
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
	for (unsigned i = 0; i < this->inputData.size(); i++)
		this->search_space.push_back(i);
	this->distMatrix = (((alphaComplex<alphaNode> *)inData.complex)->distMatrix);
	std::set<std::set<unsigned>> dsimplexes;
	std::map<std::set<unsigned>, unsigned> inner_d_1_shell;
	std::map<std::set<unsigned>, unsigned> outer_d_1_shell;
	std::set<std::set<unsigned>> simplex = {{116, 128, 167}};
	for (auto new_simplex : simplex)
	{
		dsimplexes.insert(new_simplex);
		for (auto i : new_simplex)
		{
			std::set<unsigned> key = new_simplex;
			key.erase(i);
			inner_d_1_shell[key] = i;
		}
	}
	int new_point;
	while (inner_d_1_shell.size() != 0)
	{
		for (auto iter : inner_d_1_shell)
		{
			new_point = expand_d_minus_1_simplex(iter.first, iter.second);
			if (new_point == -1)
				continue;
			std::set<unsigned> new_simplex = iter.first;
			new_simplex.insert(new_point);
			auto x = dsimplexes.insert(new_simplex);
			if (x.second == false)
				continue;
			for (auto i : iter.first)
			{
				std::cout << i << " ";
			}
			std::cout <<iter.second<< ";";
			for (auto i : new_simplex)
			{
				std::cout << i << " ";
			}
			std::cout << std::endl;
			for (auto i : new_simplex)
			{
				std::set<unsigned> key = new_simplex;
				key.erase(i);
				if (outer_d_1_shell.find(key) != outer_d_1_shell.end())
				{
					outer_d_1_shell.erase(key);
					continue;
				}
				outer_d_1_shell[key] = i;
			}
		}
		for (auto simp : inner_d_1_shell)
		{
			outer_d_1_shell.erase(simp.first);
		}

		inner_d_1_shell = outer_d_1_shell;
		outer_d_1_shell.clear();
	}
	/* 	for (auto simp : dsimplexes)
		{
			for (auto i : simp)
			{
				std::cout << i << " ";
			}
			std::cout << std::endl;
		} */
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