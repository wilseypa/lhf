#include <incrementalPipe.hpp>
#include <Eigen/Dense>
#include <limits>
#include <omp.h>
#include <random>
#include <chrono>

template <typename T>
T dot(const std::vector<T> &a, const std::vector<T> &b) { return std::inner_product(a.begin(), a.end(), b.begin(), static_cast<T>(0)); } // Dot product of two vectors

// Solution to equation of a hyperplane
template <typename nodeType>
std::vector<double> incrementalPipe<nodeType>::solvePlaneEquation(const std::vector<short> &points)
{
	int numPoints = points.size();
	Eigen::MatrixXd A(numPoints, numPoints);
	for (int i = 0; i < numPoints; i++)
		for (int j = 0; j < numPoints; j++)
			A(i, j) = inputData[points[i]][j];
	Eigen::VectorXd coefficients = A.completeOrthogonalDecomposition().solve(Eigen::VectorXd::Ones(numPoints));
	return std::vector<double>(coefficients.data(), coefficients.data() + coefficients.size());
}

// Compute the seed simplex for initialization of algorithm
template <typename nodeType>
std::vector<short> incrementalPipe<nodeType>::first_simplex()
{
	std::vector<short> simplex;
	if (this->data_set_size < this->dim + 2)
		return simplex; // Not enough points
	for (short i = 0; i < this->dim; i++)
		simplex.push_back(i); // Pseudo Random initialization of Splitting Hyperplane
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	auto rng = std::default_random_engine{};
	std::vector<short> outer_points;
	do
	{
		simplex.insert(simplex.end(), outer_points.begin(), outer_points.end());
		std::shuffle(simplex.begin(), simplex.end(), rng);
		outer_points.clear();
		simplex.erase(simplex.begin() + this->dim, simplex.end());
		auto equation = solvePlaneEquation(simplex);
		for (unsigned i = 0; i < this->data_set_size; i++)
			if (std::find(simplex.begin(), simplex.end(), i) == simplex.end() && dot(this->inputData[i], equation) > 1)
				outer_points.push_back(i);
	} while (!outer_points.empty()); // Converge Hyperplane to Convex Hull
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
void incrementalPipe<nodeType>::cospherical_handler(std::vector<short> &simp, int &tp, short &omission, std::vector<std::vector<double>> &distMatrix)
{
	auto triangulation_point = tp;
	auto temp = simp;
	temp.push_back(triangulation_point);
	std::sort(temp.begin(), temp.end());
	this->dsimplexes.insert(temp);
	for (auto i : temp)
	{
		auto new_face = temp;
		new_face.erase(std::find(new_face.begin(), new_face.end(), i));
		if (simp == new_face)
			continue;
		auto new_point = expand_d_minus_1_simplex(new_face, i, distMatrix);
		if (new_point == -1)
			continue;
		new_face.push_back(new_point);
		std::sort(new_face.begin(), new_face.end());
		this->dsimplexes.insert(new_face);
	}
}

// Compute the P_newpoint for provided Facet, P_context Pair
template <typename nodeType>
int incrementalPipe<nodeType>::expand_d_minus_1_simplex(std::vector<short> &simp, short &omission, std::vector<std::vector<double>> &distMatrix)
{
	auto normal = solvePlaneEquation(simp);
	auto p1 = utils::circumCenter(simp, this->inputData);
	bool direction = (dot(normal, this->inputData[omission]) > 1);
	std::vector<double> radius_vec(this->data_set_size, 0);
	double largest_radius = 0, smallest_radius = std::numeric_limits<double>::max() / 2, ring_radius = utils::vectors_distance(p1, inputData[simp[0]]);
	int triangulation_point = -1;
	bool flag = true;
	std::set<short> simplex(simp.begin(), simp.end());
	simplex.insert(omission);
	for (auto &new_point : this->search_space)
	{
		if (simplex.find(new_point) == simplex.end() && (direction ^ dot(normal, inputData[new_point]) > 1))
		{
			simp.push_back(new_point);
			auto temp = utils::vectors_distance(p1, inputData[new_point]);
			if (ring_radius > temp)
			{
				radius_vec[new_point] = sqrt(utils::circumRadius(simp, distMatrix));
				if (largest_radius < radius_vec[new_point])
				{
					largest_radius = radius_vec[new_point];
					triangulation_point = new_point;
					flag = false;
				}
			}
			else if (flag && 2 * smallest_radius > temp) // reduce search space
			{
				radius_vec[new_point] = sqrt(utils::circumRadius(simp, distMatrix));
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
			// cospherical_handler(simp,triangulation_point,omission,distMatrix);
			return -1;
		}
	}
	return triangulation_point;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void incrementalPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	this->inputData = inData.inputData;
	this->dim = inputData[0].size();
	this->data_set_size = inputData.size();
	for (unsigned i = 0; i < inputData.size(); i++)
		this->search_space.push_back(i);
	this->dsimplexes = {first_simplex()};
	for (auto &new_simplex : dsimplexes)
	{
		for (auto &i : new_simplex)
		{
			std::vector<short> key = new_simplex;
			key.erase(std::find(key.begin(), key.end(), i));
			this->inner_d_1_shell.emplace(key, i);
		}
	}
	short new_point;
	while (!this->inner_d_1_shell.empty())
	{
		auto iter = this->inner_d_1_shell.begin();
		std::vector<short> first_vector = iter->first;
		short omission = iter->second;
		this->inner_d_1_shell.erase(iter);
		new_point = expand_d_minus_1_simplex(first_vector, omission, inData.distMatrix);
		if (new_point == -1)
			continue;
		first_vector.push_back(new_point);
		std::sort(first_vector.begin(), first_vector.end());
		if (!this->dsimplexes.insert(first_vector).second)
			continue;
		for (auto &i : first_vector)
		{
			if (i == new_point)
				continue;
			std::vector<short> key = first_vector;
			key.erase(std::find(key.begin(), key.end(), i));
			auto temp = this->inner_d_1_shell.emplace(key, i);
			if (!temp.second)
				this->inner_d_1_shell.erase(temp.first); // Create new shell and remove collided faces max only 2 can occur.
		}
	}
	std::cout << dsimplexes.size() << std::endl;
	return;
}

// basePipe constructor
template <typename nodeType>
incrementalPipe<nodeType>::incrementalPipe() : inputData(*(new std::vector<std::vector<double>>()))
{
	this->pipeType = "incrementalPipe";
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