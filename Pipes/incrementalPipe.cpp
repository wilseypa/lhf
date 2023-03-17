#include "incrementalPipe.hpp"
#include "alphaComplex.hpp"
#include "utils.hpp"
#include <fstream>
#include <omp.h>
#include <unistd.h>

template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> c(a.size());
    for (int i = 0; i < a.size(); ++i)
        c[i] = a[i] - b[i];
    return c;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    std::cout << "[";
    for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
        std::cout << *it;
        if (it != v.end() - 1) 
            std::cout << ", ";
    }
    std::cout << "]";
    return os;
}

template <typename T>
std::set<T> operator+(std::set<T> a, T b) {
	a.insert(b);
    return a;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double c=0;
    for (int i = 0; i < a.size(); ++i)
        c += a[i] * b[i];
    return c;
}
template <typename nodeType>
int incrementalPipe<nodeType>::expand_d_minus_1_simplex(std::set<unsigned> simp, unsigned omission)
{
	auto p1=utils::circumCenter(simp, this->inputData);
	auto normal=utils::circumCenter(simp+omission, this->inputData)-p1;
	auto direction=dot(normal,this->inputData[omission]-p1)>0;
	double rad_vector[this->inputData.size()];
	double smallest_radius = 1e+308, ring_radius;
	int triangulation_point = -1;
	for (auto new_point : this->search_space)
	{
		if (simp.find(new_point) != simp.end() || new_point == omission || (direction^(dot(normal,this->inputData[new_point]-p1)<0)))
			continue;
		rad_vector[new_point] = sqrt(utils::circumRadius(simp+new_point,this->distMatrix));
		if (smallest_radius > rad_vector[new_point])
		{
			smallest_radius = rad_vector[new_point];
			triangulation_point = new_point;
		}
	}
	if (triangulation_point != -1)
	{
		auto smallest_circle_center = utils::circumCenter(simp+(unsigned)triangulation_point, this->inputData);
		ring_radius = smallest_radius;
		for (const auto& new_point : this->search_space)
		{
			if (simp.find(new_point) == simp.end() && new_point != omission && smallest_radius < rad_vector[new_point] && ring_radius > utils::vectors_distance(smallest_circle_center, this->inputData[new_point]))
			{
				smallest_radius = rad_vector[new_point];
				triangulation_point = new_point;
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
	std::set<std::set<unsigned>> simplex = {{2, 37, 73, 77, 168, 190}};
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
			auto new_simplex = iter.first;
			new_simplex=new_simplex+(unsigned)new_point;
			auto x = dsimplexes.insert(new_simplex);
			if (x.second == false)
				continue;
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
			outer_d_1_shell.erase(simp.first);
		inner_d_1_shell = outer_d_1_shell;
		outer_d_1_shell.clear();
	}
	std::cout<<dsimplexes.size();
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