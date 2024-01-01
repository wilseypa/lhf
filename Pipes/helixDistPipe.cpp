#include <helixDistPipe.hpp>
#include "utils.hpp"
#include <multifileops.hpp>
#include <fileops.hpp>

#define WRITE_CSV_OUTPUTS

/* int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <input_filename>" << std::endl;
		return 1;
	}
	MPI_Init(&argc, &argv);
	auto start_time = std::chrono::high_resolution_clock::now();
	helixDistPipe DT(argv[1]); // Use argv[1] as the input filename
	DT.run();
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time = end_time - start_time;
	MPI_Finalize();
	std::cout << "Elapsed time: " << elapsed_time.count() << " seconds" << std::endl;
	return 0;
} */

template <typename T>
static inline T dot(const std::vector<T> &a, const std::vector<T> &b) { return std::inner_product(a.begin(), a.end(), b.begin(), static_cast<T>(0)); } // Dot product of two vectors

template <typename nodeType>
helixDistPipe<nodeType>::helixDistPipe() : dim(0), data_set_size(0)
{
	this->pipeType = "delaunayPipe";
	return;
}

template <typename nodeType>
void helixDistPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	inputData = &inData.inputData;
	distMatrix = &inData.distMatrix;
	dim = inputData.size() > 0 ? inputData[0].size() : 0;
	data_set_size = inputData.size();
	for (unsigned i = 0; i < inputData.size(); i++)
		this->search_space.insert(i);
	run();
}
// Solution to equation of a hyperplane
template <typename nodeType>
inline std::vector<double> helixDistPipe<nodeType>::solvePlaneEquation(const std::vector<short> &points)
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
std::vector<short> helixDistPipe<nodeType>::first_simplex()
{
	std::vector<short> simplex;
	if (data_set_size < dim + 2)
		return simplex; // Not enough points
	for (short i = 0; i < dim; i++)
		simplex.push_back(i); // Pseudo Random initialization of Splitting Hyperplane
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	auto rng = std::default_random_engine{};
	std::vector<short> outer_points;
	do
	{
		simplex.insert(simplex.end(), outer_points.begin(), outer_points.end());
		std::shuffle(simplex.begin(), simplex.end(), rng);
		outer_points.clear();
		simplex.erase(simplex.begin() + dim, simplex.end());
		auto equation = solvePlaneEquation(simplex);
		for (unsigned i = 0; i < data_set_size; i++)
			if (std::find(simplex.begin(), simplex.end(), i) == simplex.end() && dot(inputData[i], equation) > 1)
				outer_points.push_back(i);
	} while (!outer_points.empty()); // Converge Hyperplane to Convex Hull
	double radius = 0;
	std::vector<double> center;
	for (short i = 0; i < data_set_size; i++) // BruteForce to Find last point for construction of simplex.
	{
		if (std::find(simplex.begin(), simplex.end(), i) != simplex.end())
			continue;
		simplex.push_back(i);
		center = utils::circumCenter(simplex, inputData);
		radius = utils::vectors_distance(center, inputData[i]);
		short point;
		for (point = 0; point < data_set_size; point++)
		{
			if (std::find(simplex.begin(), simplex.end(), point) == simplex.end() && utils::vectors_distance(center, inputData[point]) < radius)
				break;
		}
		if (point == data_set_size)
			break;
		simplex.pop_back();
	}
	std::sort(simplex.begin(), simplex.end());
	return simplex;
}

// Perform DFS walk on the cospherical region
template <typename nodeType>
void helixDistPipe<nodeType>::cospherical_handler(std::vector<short> &simp, int &tp, short &omission)
{
	auto triangulation_point = tp;
	auto temp = simp;
	temp.push_back(triangulation_point);
	std::sort(temp.begin(), temp.end());
	spherical_dsimplexes.insert(temp);
	for (auto i : temp)
	{
		auto new_face = temp;
		new_face.erase(std::find(new_face.begin(), new_face.end(), i));
		if (simp == new_face)
			continue;
		auto new_point = expand_d_minus_1_simplex(new_face, i);
		if (new_point == -1)
			continue;
		new_face.push_back(new_point);
		std::sort(new_face.begin(), new_face.end());
		spherical_dsimplexes.insert(new_face);
	}
}

#if 1
// Compute the P_newpoint for provided Facet, P_context Pair
template <typename nodeType>
short helixDistPipe<nodeType>::expand_d_minus_1_simplex(std::vector<short> &simp, short &omission)
{
	auto normal = solvePlaneEquation(simp);
	auto p1 = utils::circumCenter(simp, inputData);
	bool direction = (dot(normal, inputData[omission]) > 1);
	std::vector<double> radius_vec(data_set_size, 0);
	double largest_radius = 0, smallest_radius = std::numeric_limits<double>::max() / 2, ring_radius = utils::vectors_distance(p1, inputData[simp[0]]);
	int triangulation_point = -1;
	bool flag = true;
	std::set<short> simplex(simp.begin(), simp.end());
	simplex.insert(omission);
	for (auto &new_point : search_space)
	{
		if (simplex.find(new_point) == simplex.end() && (direction ^ dot(normal, inputData[new_point]) > 1))
		{
			simp.push_back(new_point);
			auto temp = utils::vectors_distance(p1, inputData[new_point]);
			if (ring_radius > temp)
			{
				auto temp_radius = utils::circumRadius(simp, distMatrix);
				radius_vec[new_point] = (temp_radius >= 0) ? sqrt(temp_radius) : utils::vectors_distance(utils::circumCenter(simp, inputData), inputData[new_point]);
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
				radius_vec[new_point] = (temp_radius >= 0) ? sqrt(temp_radius) : utils::vectors_distance(utils::circumCenter(simp, inputData), inputData[new_point]);
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

#else
// Compute the P_newpoint for provided Facet, P_context Pair
short helixDistPipe::expand_d_minus_1_simplex(std::vector<short> &simp, short &omission)
{
	// Default values
	short triangulation_point = -1;
	double largest_radius = 0;
	double smallest_radius = std::numeric_limits<double>::max() / 2;

	// Modify search_spaceF
	search_space.erase(omission);
	for (const auto &point : simp)
		search_space.erase(point);

	// Mandatory requirements
	auto normal = solvePlaneEquation(simp);
	bool direction = dot(normal, inputData[omission]) > 1;

	auto center = circumCenter(simp, inputData);

	double ring_radius = vectors_distance(center, inputData[simp[0]]);

	bool flag = true;
	return triangulation_point;
}
#endif

void reduce(std::set<std::vector<short>> &outer_dsimplexes, std::map<std::vector<short>, short> &inner_d_1_shell, std::vector<std::vector<short>> &dsimplexes)
{
	std::map<std::vector<short>, short> outer_d_1_shell;
	for (auto &new_simplex : outer_dsimplexes)
	{
		dsimplexes.push_back(new_simplex);
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
	std::for_each(std::execution::par_unseq, inner_d_1_shell.begin(), inner_d_1_shell.end(), [&](const auto &simp)
				  { outer_d_1_shell.erase(simp.first); }); // Remove faces from previous iteration
	inner_d_1_shell.clear();
	std::swap(inner_d_1_shell, outer_d_1_shell);
	return;
}

// run -> Run the configured functions of this line segment
template <typename nodeType>
void helixDistPipe<nodeType>::run()
{
	int numProcesses, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Initialization of algorithm by serial processing
	if (rank == 0)
	{
		// Refresh the directories
		std::clog << "MPI configured with " << numProcesses << " processes" << std::endl;
		std::filesystem::exists("output") && std::filesystem::remove_all("output");
		std::filesystem::exists("intermediate") && std::filesystem::remove_all("intermediate");
		std::filesystem::exists("input") && std::filesystem::remove_all("input");
		std::filesystem::create_directory("output") && std::filesystem::create_directory("intermediate") && std::filesystem::create_directory("input");
		// Perform initial iteration normally
		std::vector<std::vector<short>> initial_dsimplexes = {first_simplex()};
		std::map<std::vector<short>, short> inner_d_1_shell;
		for (auto &new_simplex : initial_dsimplexes)
		{
			for (auto &i : new_simplex)
			{
				std::vector<short> key = new_simplex;
				key.erase(std::find(key.begin(), key.end(), i));
				inner_d_1_shell.emplace(key, i);
			}
		}
		// Compute 10000 facets per process to reduce no of iterations
		for (int i = 0; i < 3; i++)
		{
			std::set<std::vector<short>> outer_dsimplexes;
			for (auto &[facet, point] : inner_d_1_shell)
			{
				std::vector<short> first_vector = facet;
				short new_point = expand_d_minus_1_simplex(first_vector, point);
				if (new_point == -1)
					continue;
				first_vector.push_back(new_point);
				std::sort(first_vector.begin(), first_vector.end());
				outer_dsimplexes.insert(first_vector);
			}
			reduce(outer_dsimplexes, inner_d_1_shell, initial_dsimplexes);
		}
		std::clog << "Iter 0 on process " << rank << " Found " << initial_dsimplexes.size() << " dsimplexes" << std::endl;
#ifdef WRITE_CSV_OUTPUTS
		std::sort(initial_dsimplexes.begin(), initial_dsimplexes.end());
		writeBinaryFile(initial_dsimplexes, "output/0_0.bin");
#endif
		initial_dsimplexes.clear();
		if (inner_d_1_shell.empty())
			return;
		writeBinaryFile(inner_d_1_shell, "input/1.dat");
	}

	// Restrict other processes from proceeding until serial task is completed
	MPI_Barrier(MPI_COMM_WORLD);

	// Compute layer by layer next simplexes
	int iter_counter = 1;
	while (true)
	{
		std::map<std::vector<short>, short> local_d_1_shell_map = readBinaryMapFile("input/" + std::to_string(iter_counter) + ".dat", numProcesses, rank);
		if (local_d_1_shell_map.empty())
			break;
		std::vector<std::vector<short>> local_dsimplexes_output;
		while (!local_d_1_shell_map.empty()) // Compute dsimplexes for individual process
		{
			auto iter = local_d_1_shell_map.begin();
			std::vector<short> first_vector = iter->first;
			short omission = iter->second;
			local_d_1_shell_map.erase(iter);
			short new_point = expand_d_minus_1_simplex(first_vector, omission);
			if (new_point == -1)
				continue;
			first_vector.insert(std::lower_bound(first_vector.begin(), first_vector.end(), new_point), new_point);
			for (auto &i : first_vector)
			{
				std::vector<short> key = first_vector;
				key.erase(std::find(key.begin(), key.end(), i));
				local_d_1_shell_map.erase(key);
			}
			local_dsimplexes_output.push_back(std::move(first_vector));
		}

		std::clog << "Iter " << iter_counter << " on process " << rank << " Found " << local_dsimplexes_output.size() << " dsimplexes" << std::endl;
		// Commit dsimplexes to file

#ifdef WRITE_CSV_OUTPUTS
		if (!local_dsimplexes_output.empty())
		{
			std::sort(local_dsimplexes_output.begin(), local_dsimplexes_output.end());
			writeBinaryFile(local_dsimplexes_output, "output/" + std::to_string(iter_counter) + "_" + std::to_string(rank) + ".bin");
		}
#endif
		// Compute facets from current layer and store to intermediate file
		std::map<std::vector<short>, short> outer_d_1_shell;
		for (auto &new_simplex : local_dsimplexes_output)
		{
			for (short i = 0; i < new_simplex.size(); i++)
			{
				std::vector<short> key = new_simplex;
				key.erase(key.begin() + i);
				auto it = outer_d_1_shell.try_emplace(std::move(key), new_simplex[i]);
				it.first->second = (it.second ? new_simplex[i] : -1);
			}
		}
		local_dsimplexes_output.clear();
		writeBinaryFile(outer_d_1_shell, "intermediate/" + std::to_string(rank) + ".dat");
		outer_d_1_shell.clear();

		// Wait for all process to commit facets
		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0)
		{
			// Perform custom multifile sort on the intermediate simplexes also remove duplicate entries
			MultiFIle<MapBinaryFile, std::pair<std::vector<short>, short>> facets("intermediate");
			facets.compressMap("input/" + std::to_string(iter_counter + 1) + ".dat", iter_counter);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		std::filesystem::remove("intermediate/" + std::to_string(rank) + ".dat");
		iter_counter++;
	}
#ifdef WRITE_CSV_OUTPUTS
	if (rank == 0)
	{
		MultiFIle<VectorBinaryFile, std::vector<short>> dsimplexes("output");
		std::clog << "Found " << dsimplexes.writeCSV("dsimplexes.csv") << " dsimplexes for datset" << std::endl;
	}
#endif

	return;
}

template class helixDistPipe<simplexNode>;
template class helixDistPipe<alphaNode>;
template class helixDistPipe<witnessNode>;