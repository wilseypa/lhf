#include "helixPipe.hpp"
#include <Eigen/Dense>
#include <limits>
#include <omp.h>
#include <random>
#include <chrono>
#include <execution>
#include <fstream>
#include <ranges>
// #define PARALLEL

// Solution to equation of a hyperplane
template <typename nodeType>
std::vector<double> helixPipe<nodeType>::solvePlaneEquation(const std::vector<short> &points)
{
	int numPoints = points.size();
	Eigen::MatrixXd A(numPoints, numPoints);
	for (int i = 0; i < numPoints; i++)
		A.row(i) = Eigen::Map<Eigen::VectorXd>(this->inputData[points[i]].data(), numPoints);
	Eigen::VectorXd coefficients = A.completeOrthogonalDecomposition().solve(Eigen::VectorXd::Ones(numPoints));
	return std::vector<double>(coefficients.data(), coefficients.data() + coefficients.size());
}

// Compute the seed simplex for initialization of algorithm
template <typename nodeType>
std::vector<short> helixPipe<nodeType>::first_simplex()
{
	std::vector<short> simplex(this->dim);
	std::iota(simplex.begin(), simplex.end(), 0); // Pseudo Random initialization of Splitting Hyperplane
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	auto rng = std::default_random_engine{};
	std::vector<short> outer_points;
	while (true)
	{
		simplex.insert(simplex.end(), outer_points.begin(), outer_points.end());
		std::shuffle(simplex.begin(), simplex.end(), rng);
		outer_points.clear();
		simplex.resize(this->dim);
		auto equation = this->solvePlaneEquation(simplex);
		for (unsigned i = 0; i < this->data_set_size; i++)
			if (std::find(simplex.begin(), simplex.end(), i) == simplex.end() && utils::dot(this->inputData[i], equation) > 1)
				outer_points.push_back(i);
		if (outer_points.empty()) // Converge Hyperplane to Convex Hull
			break;
	}
	double radius = 0;
	std::vector<double> center;
	for (short i = 0; i < this->data_set_size; i++) // BruteForce to Find last point for construction of simplex.
	{
		if (std::find(simplex.begin(), simplex.end(), i) != simplex.end())
			continue;
		simplex.push_back(i);
		center = utils::circumCenter(simplex, this->inputData);
		radius = utils::vectors_distance(center, this->inputData[i]);
		short point;
		for (point = 0; point < this->data_set_size; point++)
		{
			if (std::find(simplex.begin(), simplex.end(), point) == simplex.end() && utils::vectors_distance(center, this->inputData[point]) < radius)
				break;
		}
		if (point == this->data_set_size)
			break;
		simplex.pop_back();
	}
	std::sort(simplex.begin(), simplex.end());
	return simplex;
}

// Perform DFS walk on the cospherical region
template <typename nodeType>
void helixPipe<nodeType>::cospherical_handler(std::vector<short> &simp, int &tp, short &omission, std::vector<std::vector<double>> &distMatrix)
{
	auto triangulation_point = tp;
	auto temp = simp;
	temp.push_back(triangulation_point);
	std::sort(temp.begin(), temp.end());
	this->spherical_dsimplexes.insert(temp);
	for (auto i : temp)
	{
		auto new_face = temp;
		new_face.erase(std::find(new_face.begin(), new_face.end(), i));
		if (simp == new_face)
			continue;
		auto new_point = this->expand_d_minus_1_simplex(new_face, i, distMatrix);
		if (new_point == -1)
			continue;
		new_face.push_back(new_point);
		std::sort(new_face.begin(), new_face.end());
		this->spherical_dsimplexes.insert(new_face);
	}
}

#if 0
// Compute the P_newpoint for provided Facet, P_context Pair
template <typename nodeType>
short helixPipe<nodeType>::expand_d_minus_1_simplex(std::vector<short> &simp, short &omission, std::vector<std::vector<double>> &distMatrix)
{
	auto normal = solvePlaneEquation(simp);
	auto p1 = utils::circumCenter(simp, this->inputData);
	bool direction = (dot(normal, this->inputData[omission]) > 1);
	std::vector<double> radius_vec(this->data_set_size, 0);
	double largest_radius = 0, smallest_radius = std::numeric_limits<double>::max() / 2, ring_radius = utils::vectors_distance(p1, this->inputData[simp[0]]);
	int triangulation_point = -1;
	bool flag = true;
	std::set<short> simplex(simp.begin(), simp.end());
	simplex.insert(omission);
	for (auto &new_point : this->search_space)
	{
		if (simplex.find(new_point) == simplex.end() && (direction ^ dot(normal, this->inputData[new_point]) > 1))
		{
			simp.push_back(new_point);
			auto temp = utils::vectors_distance(p1, this->inputData[new_point]);
			if (ring_radius > temp)
			{
				auto temp_radius = utils::circumRadius(simp, distMatrix);
				radius_vec[new_point] = (temp_radius >= 0) ? sqrt(temp_radius) : utils::vectors_distance(utils::circumCenter(simp, this->inputData), this->inputData[new_point]);
				if (largest_radius < radius_vec[new_point])
				{
					largest_radius = radius_vec[new_point];
					triangulation_point = new_point;
					flag = false;
				}
			}
			else if (flag && 2 * smallest_radius > temp) // reduce search space
			{
				auto temp_radius = utils::circumRadius(simp, distMatrix);
				radius_vec[new_point] = (temp_radius >= 0) ? sqrt(temp_radius) : utils::vectors_distance(utils::circumCenter(simp, this->inputData), this->inputData[new_point]);
				if (smallest_radius > radius_vec[new_point])
				{
					smallest_radius = radius_vec[new_point];
					triangulation_point = new_point;
				}
			}
			simp.pop_back();
		}
	}
	int count = 0;
	if (triangulation_point != -1)
	{
		double triangulation_radius = radius_vec[triangulation_point];
		count = std::count_if(radius_vec.begin(), radius_vec.end(), [triangulation_radius](double val)
							  { return std::abs(1 - (val / triangulation_radius)) <= 0.000000000001; });
		if (count != 1)
		{
			std::cout << "Triangulation radius is " << triangulation_radius << std::endl;
#ifdef PARALLEL
#pragma omp critical
			{
				// cospherical_handler(simp,triangulation_point,omission,distMatrix);
			}
#else
			// cospherical_handler(simp,triangulation_point,omission,distMatrix);
#endif
			return -1;
		}
	}
	return triangulation_point;
}

#else
// Compute the P_newpoint for provided Facet, P_context Pair
template <typename nodeType>
short helixPipe<nodeType>::expand_d_minus_1_simplex(std::vector<short> &simp, short &omission, std::vector<std::vector<double>> &distMatrix)
{
	// Mandatory requirements
	auto normal = this->solvePlaneEquation(simp);
	bool direction = utils::dot(normal, inputData[omission]) > 1;

	// Filter points opposite to hyperplane
	for (size_t idx = 0; idx < this->active_data.size(); ++idx)
		this->active_data[idx] = direction ^ (utils::dot(normal, inputData[idx]) > 1) ? 0 : std::numeric_limits<double>::infinity();

	// Deactivate points from original simplex
	this->active_data[omission] = std::numeric_limits<double>::infinity();
	for (const auto &point : simp)
		this->active_data[point] = std::numeric_limits<double>::infinity();

	// Early return if plane is on convex hull
	if (std::all_of(this->active_data.begin(), this->active_data.end(), [](double isActive)
					{ return isActive == std::numeric_limits<double>::infinity(); }))
		return -1;

	// Default values
	short triangulation_point = -1;
	// Calculate the circumcenter of the facet
	auto center = utils::circumCenter(simp, inputData);

	for (size_t idx = 0; idx < this->active_data.size(); ++idx)
		if (this->active_data[idx] != std::numeric_limits<double>::infinity())
			this->active_data[idx] = utils::vectors_distance(center, inputData[idx]);

	auto ring_radius = utils::vectors_distance(center, this->inputData[simp[0]]);

	// Check for any distances less than the distance to the first simplex point
	if (std::find_if(this->active_data.begin(), this->active_data.end(), [&](double distance)
					 { return distance < ring_radius; }) == this->active_data.end())
	{
		// Search outside space with smallest circumradius
		double smallest_radius = std::numeric_limits<double>::max() / 2;
		for (size_t idx = 0; idx < this->active_data.size(); ++idx)
		{
			if (this->active_data[idx] != std::numeric_limits<double>::infinity() && this->active_data[idx] > ring_radius && this->active_data[idx] < 2 * smallest_radius)
			{
				simp.push_back(idx);
				auto temp_radius = utils::circumRadius(simp, distMatrix);
				this->active_data[idx] = (temp_radius >= 0) ? sqrt(temp_radius) : utils::vectors_distance(utils::circumCenter(simp, this->inputData), this->inputData[idx]);
				simp.pop_back();
				if (smallest_radius > this->active_data[idx])
				{
					smallest_radius = this->active_data[idx];
					triangulation_point = idx;
				}
			}
			else
				this->active_data[idx] = 0;
		}
	}
	else
	{
		// Search inside space with largest circumradius
		double largest_radius = 0;
		for (size_t idx = 0; idx < this->active_data.size(); ++idx)
		{
			if (this->active_data[idx] != std::numeric_limits<double>::infinity() && this->active_data[idx] < ring_radius)
			{
				simp.push_back(idx);
				auto temp_radius = utils::circumRadius(simp, distMatrix);
				this->active_data[idx] = (temp_radius >= 0) ? sqrt(temp_radius) : utils::vectors_distance(utils::circumCenter(simp, this->inputData), this->inputData[idx]);
				simp.pop_back();
				if (largest_radius < this->active_data[idx])
				{
					largest_radius = this->active_data[idx];
					triangulation_point = idx;
				}
			}
			else
				this->active_data[idx] = 0;
		}
	}
	double triangulation_radius = this->active_data[triangulation_point];
	int count = std::count_if(this->active_data.begin(), this->active_data.end(), [triangulation_radius](double val)
							  { return std::abs(1 - (val / triangulation_radius)) <= 0.000000000001; });
	if (count != 1)
	{
		std::cout << "Cospherical Region Found at triangulation radius " << triangulation_radius << std::endl;
#ifdef PARALLEL
#pragma omp critical
		{
			// this->cospherical_handler(simp,triangulation_point,omission,distMatrix);
		}
#else
		// this->cospherical_handler(simp,triangulation_point,omission,distMatrix);
#endif
		return -1;
	}
	return triangulation_point;
}
#endif

template <typename nodeType>
void helixPipe<nodeType>::reduce(std::set<std::vector<short>> &outer_dsimplexes, std::vector<std::pair<std::vector<short>, short>> &inner_d_1_shell, std::vector<std::vector<short>> &dsimplexes)
{
	std::map<std::vector<short>, short> outer_d_1_shell;
	for (auto &new_simplex : outer_dsimplexes)
	{
		dsimplexes.emplace_back(new_simplex);
		for (short i = 0; i < new_simplex.size(); i++)
		{
			std::vector<short> key = new_simplex;
			key.erase(key.begin() + i);
			auto it = outer_d_1_shell.try_emplace(std::move(key), new_simplex[i]);
			if (!it.second)
				outer_d_1_shell.erase(it.first); // Create new shell and remove collided faces max only 2 can occur.
		}
	}
	outer_dsimplexes.clear();
#ifndef NO_PARALLEL_ALGORITHMS
	std::for_each(std::execution::par_unseq, inner_d_1_shell.begin(), inner_d_1_shell.end(), [&](const auto &simp)
				  { outer_d_1_shell.erase(simp.first); }); // Remove faces from previous iteration
#else
	std::for_each(inner_d_1_shell.begin(), inner_d_1_shell.end(), [&](const auto &simp)
				  { outer_d_1_shell.erase(simp.first); });
#endif
	inner_d_1_shell.clear();
	inner_d_1_shell.reserve(outer_d_1_shell.size());
	std::move(outer_d_1_shell.begin(), outer_d_1_shell.end(), std::back_inserter(inner_d_1_shell));
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void helixPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	this->inputData = inData.inputData;
	this->data_set_size = this->inputData.size();
	if (!this->data_set_size)
		return; // Empty inputData
	this->dim = this->inputData.size() > 0 ? this->inputData[0].size() : 0;
	if (this->data_set_size < this->dim + 2)
		return; // Not enough points

	this->active_data.resize(this->data_set_size, 0);

	this->dsimplexmesh = {this->first_simplex()};
	short new_point;
#ifdef PARALLEL
	std::vector<std::pair<std::vector<short>, short>> inner_d_1_shell;
	for (auto &new_simplex : this->dsimplexmesh)
	{
		for (auto &i : new_simplex)
		{
			std::vector<short> key = new_simplex;
			key.erase(std::find(key.begin(), key.end(), i));
			inner_d_1_shell.push_back(std::make_pair(key, i));
		}
	}
	std::set<std::vector<short>> outer_dsimplexes;
	while (inner_d_1_shell.size() != 0)
	{
#pragma omp parallel for private(new_point, this->this->active_data)
		for (int i = 0; i < inner_d_1_shell.size(); i++)
		{
			auto iter = inner_d_1_shell[i];
			new_point = this->expand_d_minus_1_simplex(iter.first, iter.second, inData.distMatrix);
			if (new_point == -1)
				continue;
			iter.first.push_back(new_point);
			std::sort(iter.first.begin(), iter.first.end());
			if (outer_dsimplexes.find(iter.first) == outer_dsimplexes.end())
#pragma omp critical
				outer_dsimplexes.insert(iter.first);
		}
		std::cout << "Intermediate dsimplex size " << outer_dsimplexes.size() << std::endl;
		this->reduce(outer_dsimplexes, inner_d_1_shell, this->dsimplexmesh);
	}
#else
	std::map<std::vector<short>, short> inner_d_1_shell;
	for (auto &new_simplex : this->dsimplexmesh)
	{
		for (auto &i : new_simplex)
		{
			std::vector<short> key = new_simplex;
			key.erase(std::find(key.begin(), key.end(), i));
			inner_d_1_shell.emplace(key, i);
		}
	}
	while (!inner_d_1_shell.empty())
	{
		auto iter = inner_d_1_shell.begin();
		std::vector<short> first_vector = iter->first;
		short omission = iter->second;
		inner_d_1_shell.erase(iter);
		new_point = this->expand_d_minus_1_simplex(first_vector, omission, inData.distMatrix);
		if (new_point == -1)
			continue;
		first_vector.push_back(new_point);
		std::sort(first_vector.begin(), first_vector.end());
		this->dsimplexmesh.push_back(first_vector);
		for (auto &i : first_vector)
		{
			if (i == new_point)
				continue;
			std::vector<short> key = first_vector;
			key.erase(std::find(key.begin(), key.end(), i));
			auto temp = inner_d_1_shell.emplace(key, i);
			if (!temp.second)
				inner_d_1_shell.erase(temp.first); // Create new shell and remove collided faces max only 2 can occur.
		}
	}
#endif
	std::vector<std::vector<int>> cache(this->data_set_size, std::vector<int>(this->data_set_size, 0));
	std::sort(this->dsimplexmesh.begin(),this->dsimplexmesh.end());
	for (auto simplex : this->dsimplexmesh)
		for (int i = 0; i < simplex.size() - 1; i++)
			for (int j = i + 1; j < simplex.size(); j++)
			{
				cache[simplex[i]][simplex[j]] = 1;
				cache[simplex[j]][simplex[i]] = 1;
			}
	// Code to print the cache matrix
	for (int i = 0; i < this->data_set_size; i++)
	{
		for (int j = 0; j < this->data_set_size; j++)
		{
			std::cout << cache[i][j] << " "; // Print each element followed by a space
		}
		std::cout << std::endl; // After each row, print a newline
	}
	return;
}

// basePipe constructor
template <typename nodeType>
helixPipe<nodeType>::helixPipe()
{
	this->pipeType = "helixPipe";
	return;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool helixPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
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
	this->ut.writeDebug("helixPipe", "Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void helixPipe<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	std::ofstream file;
	file.open("output/" + this->pipeType + "_output.csv");
	// code to print the data
	for (const auto &simplex : this->dsimplexmesh)
	{
		if (simplex.empty())
			continue;
		for (size_t i = 0; i < simplex.size() - 1; i++)
			file << simplex[i] << ",";
		file << simplex.back() << "\n";
	}

	file.close();
	return;
}

template class helixPipe<simplexNode>;
template class helixPipe<alphaNode>;
template class helixPipe<witnessNode>;