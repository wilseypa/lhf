#include <incrementalPipe.hpp>
#include <Eigen/Dense>
#include <limits>
#include <omp.h>
#include <random>

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
{
	os << "[ ";
	for (const auto &elem : vec)
	{
		os << elem << ", ";
	}
	os << "]";
	return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::set<T> &set)
{
	os << "{ ";
	for (const auto &elem : set)
	{
		os << elem << ",";
	}
	os << "}";
	return os;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
	std::vector<T> temp;
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp), [](double e1, double e2)
				   { return (e1 - e2); });
	return temp;
}

template <typename T>
T dot(const std::vector<T> &a, const std::vector<T> &b)
{
	std::vector<T> temp;
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp), [](T e1, T e2)
				   { return (e1 * e2); });
	return std::accumulate(temp.begin(), temp.end(), 0.0);
}

int bruteforce(std::set<unsigned> simplex, std::vector<std::vector<double>> &inputData, short &omission)
{
	for (unsigned i = 0; i < inputData.size(); i++)
	{
		if (simplex.find(i) != simplex.end() || i == omission)
			continue;
		simplex.insert(i);
		auto center = utils::circumCenter(simplex, inputData);
		auto radius = utils::vectors_distance(center, inputData[i]);
		for (unsigned point = 0; point < inputData.size(); point++)
		{
			if (simplex.find(point) != simplex.end())
				continue;
			if (utils::vectors_distance(center, inputData[point]) < radius)
				break;
			else if (point == inputData.size() - 1)
				return i;
		}
		simplex.erase(i);
	}
	return -1;
}

std::vector<double> solvePlaneEquation(const std::vector<short> &points, const std::vector<std::vector<double>> &inputData)
{
	int numPoints = points.size();
	int numDimensions = inputData[0].size();
	Eigen::MatrixXd A(numPoints, numDimensions + 1);
	Eigen::VectorXd B(numPoints);

	for (int i = 0; i < numPoints; i++)
	{
		for (int j = 0; j < numDimensions; j++)
		{
			A(i, j) = inputData[points[i]][j];
		}
		A(i, numDimensions) = 0;
		B(i) = 1;
	}
	Eigen::VectorXd coefficients = A.completeOrthogonalDecomposition().solve(B);
	return std::vector<double>(coefficients.data(), coefficients.data() + coefficients.size());
}

int validate(std::vector<short> &simp, std::vector<std::vector<double>> &inputData, int &triangulation_point, short &omission)
{
	std::set<unsigned> simplex(simp.begin(), simp.end());
	std::vector<double> center = utils::circumCenter(simplex, inputData);
	double radius = utils::vectors_distance(center, inputData[*simp.begin()]);
	unsigned point;
	int temp;
	for (point = 0; point < inputData.size(); point++)
	{
		if (simplex.find(point) == simplex.end() && utils::vectors_distance(center, inputData[point]) < radius)
		{
			simplex.erase(triangulation_point);
			temp = bruteforce(simplex, inputData, omission);
			std::cout << temp << std::endl;
			if (temp == -1)
				return 0;
			simp.erase(std::find(simp.begin(), simp.end(), triangulation_point));
			simp.push_back(temp);
			std::sort(simp.begin(), simp.end());
			return 1;
		}
	}
	return (point == inputData.size()) ? 1 : 0;
}

std::vector<short> first_simplex(std::vector<std::vector<double>> &inputData, std::vector<std::vector<double>> &distMatrix)
{
	unsigned dim = inputData[0].size();
	unsigned n_pts = inputData.size();
	std::vector<short> simplex;
	if (n_pts < dim + 2)
		return simplex;
	while (simplex.size() != dim)
	{
		double curr_min_x = std::numeric_limits<double>::max();
		short min_x_idx = std::numeric_limits<short>::max();
		for (unsigned i = 0; i < n_pts; i++)
		{
			if (std::find(simplex.begin(), simplex.end(), i) != simplex.end())
				continue;
			double x = inputData[i][0];
			if (x < curr_min_x)
			{
				curr_min_x = x;
				min_x_idx = i;
			}
		}
		simplex.push_back(min_x_idx);
	}
	std::vector<short> outer_points;
	auto equation = solvePlaneEquation(simplex, inputData);
	for (unsigned i = 0; i < n_pts; i++)
		if (std::find(simplex.begin(), simplex.end(), i) == simplex.end() && dot(inputData[i], equation) > 1)
			outer_points.push_back(i);
	auto rng = std::default_random_engine{};
	while (outer_points.size() != 0)
	{
		for (auto i : outer_points)
			simplex.push_back(i);
		std::shuffle(simplex.begin(), simplex.end(), rng);
		outer_points.clear();
		while (simplex.size() != dim)
			simplex.pop_back();
		auto equation = solvePlaneEquation(simplex, inputData);
		for (unsigned i = 0; i < n_pts; i++)
			if (std::find(simplex.begin(), simplex.end(), i) == simplex.end() && dot(inputData[i], equation) > 1)
				outer_points.push_back(i);
	}
	double radius = 0;
	std::vector<double> center;
	std::set<unsigned> simplex_set(simplex.begin(), simplex.end());
	for (unsigned i = 0; i < inputData.size(); i++)
	{
		if (simplex_set.find(i) != simplex_set.end())
			continue;
		simplex_set.insert(i);
		center = utils::circumCenter(simplex_set, inputData);
		radius = utils::vectors_distance(center, inputData[i]);
		unsigned point;
		for (point = 0; point < inputData.size(); point++)
		{
			if (simplex_set.find(point) != simplex_set.end())
				continue;
			if (utils::vectors_distance(center, inputData[point]) < radius)
				break;
		}
		if (point == inputData.size())
		{
			simplex_set.clear();
			simplex.push_back(i);
			break;
		}
		simplex_set.erase(i);
	}
	return simplex;
}

int expand_d_minus_1_simplex(std::vector<short> &simp_vector, short &omission, std::vector<std::vector<double>> &inputData, std::vector<unsigned> &search_space, std::vector<std::vector<double>> &distMatrix)
{
	std::set<unsigned> simp(simp_vector.begin(), simp_vector.end());
	auto normal = solvePlaneEquation(simp_vector, inputData);
	auto p1 = utils::circumCenter(simp, inputData);
	auto direction = (dot(normal, inputData[omission]) > 1);
	double largest_radius = 0, curr_radius = 0, smallest_radius = std::numeric_limits<double>::max(), ring_radius = utils::vectors_distance(p1, inputData[simp_vector[0]]);
	int triangulation_point = -1;
	bool flag = true;
	for (auto &new_point : search_space)
	{
		if (simp.find(new_point) == simp.end() && new_point != omission && (direction ^ dot(normal, inputData[new_point]) > 1))
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

void reduce(std::set<std::vector<short>> &outer_dsimplexes, std::vector<std::vector<short>> &inner_d_1_shell, std::vector<std::vector<short>> &dsimplexes)
{
	std::map<std::vector<short>, short> outer_d_1_shell;
	for (auto &new_simplex : outer_dsimplexes)
	{
		// dsimplexes.push_back(new_simplex);
		for (auto &i : new_simplex)
		{
			std::vector<short> key = new_simplex;
			key.erase(std::find(key.begin(), key.end(), i));
			if (!outer_d_1_shell.emplace(key, i).second)
				outer_d_1_shell.erase(key); // Create new shell and remove collided faces max only 2 can occur.
		}
	}
	for (auto &simp : inner_d_1_shell)
	{
		simp.pop_back();
		outer_d_1_shell.erase(simp); // Remove faces from previous iteration
	}
	inner_d_1_shell.clear();
	for (auto &simp : outer_d_1_shell)
	{
		auto temp = simp.first;
		temp.push_back(simp.second);
		inner_d_1_shell.push_back(temp); // Remove faces from previous iteration
	}
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
	std::vector<std::vector<double>> distMatrix = inData.distMatrix;
	std::vector<std::vector<short>> dsimplexes = {first_simplex(inputData, distMatrix)};
	std::vector<std::vector<short>> inner_d_1_shell;
	std::set<std::vector<short>> outer_dsimplexes;
	for (auto &new_simplex : dsimplexes)
	{
		for (auto &i : new_simplex)
		{
			std::vector<short> key = new_simplex;
			key.erase(std::find(key.begin(), key.end(), i));
			key.push_back(i);
			inner_d_1_shell.push_back(key);
		}
	}
	int new_point;
	while (inner_d_1_shell.size() != 0)
	{
#pragma omp parallel for private(new_point)
		for (int i = 0; i < inner_d_1_shell.size(); i++)
		{
			auto iter = inner_d_1_shell[i];
			short omission = iter.back();
			iter.pop_back();
			new_point = expand_d_minus_1_simplex(iter, omission, inputData, search_space, distMatrix);
			if (new_point == -1)
				continue;
			iter.push_back(new_point);
			std::sort(iter.begin(), iter.end());
			if (outer_dsimplexes.find(iter) == outer_dsimplexes.end() && validate(iter, inputData, new_point, omission) == 1)
#pragma omp critical
				outer_dsimplexes.insert(iter);
		}
		std::cout << "Intermediate dsimplex size " << outer_dsimplexes.size() << std::endl;
		reduce(outer_dsimplexes, inner_d_1_shell, dsimplexes);
	}
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