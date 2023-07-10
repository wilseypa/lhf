#include <incrementalPipe.hpp>
#include <Eigen/Dense>
#include <limits>
#include <omp.h>
#include <random>
#include <chrono>

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) // Vector Subtraction
{
	std::vector<T> temp;
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp), [](double e1, double e2)
				   { return (e1 - e2); });
	return temp;
}

template <typename T>
T dot(const std::vector<T> &a, const std::vector<T> &b) //Vector Dot Product
{
	std::vector<T> temp;
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(temp), [](T e1, T e2)
				   { return (e1 * e2); });
	return std::accumulate(temp.begin(), temp.end(), 0.0);
}

int bruteforce(std::vector<short> simp, std::vector<std::vector<double>> &inputData, short &omission)
{
	std::set<unsigned> simplex(simp.begin(),simp.end());
	for (unsigned i = 0; i < inputData.size(); i++)
	{
		if (simplex.find(i) != simplex.end() || i == omission)
			continue;
		simplex.insert(i);
		auto center = utils::circumCenter(simplex, inputData);
		auto radius = utils::vectors_distance(center, inputData[i]);
		unsigned point;
		for (point = 0; point < inputData.size(); point++)
		{
			if (simplex.find(point) == simplex.end() && utils::vectors_distance(center, inputData[point]) < radius)
				break;
		}
		if (point == inputData.size())	return i;
		simplex.erase(i);
	}
	return -1;
}

template <typename nodeType>
short incrementalPipe<nodeType>::validate(std::vector<short>& simp, short &triangulation_point, short &omission)
{
	simp.push_back(triangulation_point);
	std::vector<double> center = utils::circumCenter(simp, this->inputData);
	double radius = utils::vectors_distance(center, this->inputData[*simp.begin()]);	
	for (short point = 0; point < this->data_set_size; point++)
	{
		if (std::find(simp.begin(), simp.end(), point) == simp.end() && utils::vectors_distance(center, this->inputData[point]) < 0.99999 * radius)
		{
			simp.pop_back();
			std::cout<<"Invalidated"<<std::endl;
			return bruteforce(simp, this->inputData, omission);
		}
	}
	simp.pop_back();
	return triangulation_point;
}

template <typename nodeType>
std::vector<double> incrementalPipe<nodeType>::solvePlaneEquation(const std::vector<short> &points)
{
	int numPoints = points.size();
	Eigen::MatrixXd A(numPoints, this->dim + 1);
	Eigen::VectorXd B(numPoints);
	for (int i = 0; i < numPoints; i++)
	{
		for (int j = 0; j < this->dim; j++)
			A(i, j) = this->inputData[points[i]][j];
		A(i, this->dim) = 0;
		B(i) = 1;
	}
	Eigen::VectorXd coefficients = A.completeOrthogonalDecomposition().solve(B);
	return std::vector<double>(coefficients.data(), coefficients.data() + coefficients.size());
}

template <typename nodeType>
std::vector<short> incrementalPipe<nodeType>::first_simplex()
{
	std::vector<short> simplex;
	if (this->data_set_size < this->dim + 2)
		return simplex; //Not enough points
	for(short i=0;i < this->dim;i++)
		simplex.push_back(i); //Pseudo Random initialization of Splitting Hyperplane
	std::vector<short> outer_points;
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	auto rng = std::default_random_engine{seed};
	do
	{
		for (auto i : outer_points)
			simplex.push_back(i);
		std::shuffle(simplex.begin(), simplex.end(), rng);
		outer_points.clear();
		while (simplex.size() != this->dim)
			simplex.pop_back();
		auto equation = solvePlaneEquation(simplex);
		for (unsigned i = 0; i < this->data_set_size; i++)
			if (std::find(simplex.begin(), simplex.end(), i) == simplex.end() && dot(this->inputData[i], equation) > 1)
				outer_points.push_back(i);
	} while (!outer_points.empty()); //Converge Hyperplane to Convex Hull
	double radius = 0;
	std::vector<double> center;
	std::set<unsigned> simplex_set(simplex.begin(), simplex.end());
	for (unsigned i = 0; i < this->data_set_size; i++) //BruteForce to Find last point for construction of simplex.
	{
		if (simplex_set.find(i) != simplex_set.end())
			continue;
		simplex_set.insert(i);
		center = utils::circumCenter(simplex_set, this->inputData);
		radius = utils::vectors_distance(center, this->inputData[i]);
		unsigned point;
		for (point = 0; point < this->data_set_size; point++)
		{
			if (simplex_set.find(point) == simplex_set.end() && utils::vectors_distance(center, this->inputData[point]) < radius)
				break;
		}
		if (point == this->data_set_size)
		{
			simplex_set.clear();
			simplex.push_back(i);
			break;
		}
		simplex_set.erase(i);
	}
	std::sort(simplex.begin(), simplex.end());
	return simplex;
}

template <typename nodeType>
int incrementalPipe<nodeType>::expand_d_minus_1_simplex(std::vector<short> &simp, short &omission)
{
	auto normal = solvePlaneEquation(simp);
	auto p1 = utils::circumCenter(simp, this->inputData);
	bool direction = (dot(normal, this->inputData[omission]) > 1);
	std::vector<double> radius_vec(this->data_set_size, 0);
	double largest_radius = 0, smallest_radius = std::numeric_limits<double>::max()/2, ring_radius = utils::vectors_distance(p1, inputData[simp[0]]);
	int triangulation_point = -1;
	bool flag = true;
	std::set<short> simplex(simp.begin(),simp.end());
	simplex.insert(omission);
	for (auto &new_point : this->search_space)
	{
		if (simplex.find(new_point)==simplex.end() && (direction ^ dot(normal, inputData[new_point]) > 1))
		{
			simp.push_back(new_point);
			auto temp=utils::vectors_distance(p1, inputData[new_point]);
			if (ring_radius > temp)
			{
				radius_vec[new_point] = sqrt(utils::circumRadius(simp, this->distMatrix));
				if (largest_radius < radius_vec[new_point])
				{
					largest_radius = radius_vec[new_point];
					triangulation_point = new_point;
					flag = false;
				}
			}
			else if (flag && 2*smallest_radius > temp) //reduce search space
			{
				radius_vec[new_point] = sqrt(utils::circumRadius(simp, this->distMatrix));
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
		if (count != 1){
		std::cout<< "Cospherical " << simp << " " << triangulation_point << " with "<< count << "no of points awith radius of " << triangulation_radius <<std::endl;
		return -1;
		}
	}
	return triangulation_point;
}

// basePipe constructor
template <typename nodeType>
incrementalPipe<nodeType>::incrementalPipe()  : inputData(*(new std::vector<std::vector<double>>())), distMatrix(*(new std::vector<std::vector<double>>()))
{
	this->pipeType = "incrementalPipe";
	return;
}

// runPipe -> Run the configured functions of this pipeline segment
template <typename nodeType>
void incrementalPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	this->inputData = inData.inputData;
	this->distMatrix = inData.distMatrix;
	this->dim = inputData[0].size();
	this->data_set_size = inputData.size();
	for (unsigned i = 0; i < inputData.size(); i++)
		this->search_space.push_back(i);
	std::vector<std::vector<short>> dsimplexes = {first_simplex()};
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
		new_point = expand_d_minus_1_simplex(first_vector, omission);
		if (new_point == -1)
			continue;
  		new_point = validate(first_vector, new_point, omission);
		if (new_point == -1)
			continue;  
		first_vector.push_back(new_point);
		std::sort(first_vector.begin(), first_vector.end());
		dsimplexes.push_back(first_vector);
		for (auto& i: first_vector)
		{
			if (i == new_point)
				continue;
			std::vector<short> key = first_vector;
			key.erase(std::find(key.begin(), key.end(), i));
			if (!this->inner_d_1_shell.emplace(key, i).second)
				this->inner_d_1_shell.erase(key); // Create new shell and remove collided faces max only 2 can occur.
		}
	}
	std::cout << dsimplexes.size() << std::endl;
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