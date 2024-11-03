#include "helixDistPipe.hpp"
#include "utils.hpp"
#include "readInput.hpp"
#include "writeOutput.hpp"
#include "multifileops.hpp"
#include <limits>
#include <random>
#include <chrono>
#include <execution>
#include <filesystem>
#include <Eigen/Dense>
#include <omp.h>
#include <mpi.h>

#define WRITE_CSV_OUTPUTS

template <typename nodeType>
helixDistPipe<nodeType>::helixDistPipe() : dim(0), data_set_size(0)
{
	this->pipeType = "helixDistPipe";
	return;
}

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
#ifndef NO_PARALLEL_ALGORITHMS
	std::for_each(std::execution::par_unseq, inner_d_1_shell.begin(), inner_d_1_shell.end(), [&](const auto &simp)
				  { outer_d_1_shell.erase(simp.first); }); // Remove faces from previous iteration
#else
	std::for_each(inner_d_1_shell.begin(), inner_d_1_shell.end(), [&](const auto &simp)
				  { outer_d_1_shell.erase(simp.first); });
#endif
	inner_d_1_shell.clear();
	std::swap(inner_d_1_shell, outer_d_1_shell);
	return;
}

// run -> Run the configured functions of this line segment
template <typename nodeType>
void helixDistPipe<nodeType>::runPipe(pipePacket<nodeType> &inData)
{
	inputData = inData.inputData;
	distMatrix = inData.distMatrix;
	dim = inputData.size() > 0 ? inputData[0].size() : 0;
	data_set_size = inputData.size();
	for (unsigned i = 0; i < inputData.size(); i++)
		this->search_space.insert(i);

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
		std::vector<std::vector<short>> initial_dsimplexes = {this->first_simplex()};
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
				short new_point = this->expand_d_minus_1_simplex(first_vector, point, distMatrix);
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
		writeOutput::writeBinary(initial_dsimplexes, "output/0_0.bin");
#endif
		initial_dsimplexes.clear();
		if (inner_d_1_shell.empty())
			return;
		writeOutput::writeBinary(inner_d_1_shell, "input/1.dat");
	}

	// Restrict other processes from proceeding until serial task is completed
	MPI_Barrier(MPI_COMM_WORLD);

	// Compute layer by layer next simplexes
	int iter_counter = 1;
	while (true)
	{
		std::map<std::vector<short>, short> local_d_1_shell_map = readInput::readBinaryMap("input/" + std::to_string(iter_counter) + ".dat", numProcesses, rank);
		if (local_d_1_shell_map.empty())
			break;
		std::vector<std::vector<short>> local_dsimplexes_output;
		while (!local_d_1_shell_map.empty()) // Compute dsimplexes for individual process
		{
			auto iter = local_d_1_shell_map.begin();
			std::vector<short> first_vector = iter->first;
			short omission = iter->second;
			local_d_1_shell_map.erase(iter);
			short new_point = this->expand_d_minus_1_simplex(first_vector, omission, distMatrix);
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
			writeOutput::writeBinary(local_dsimplexes_output, "output/" + std::to_string(iter_counter) + "_" + std::to_string(rank) + ".bin");
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
		writeOutput::writeBinary(outer_d_1_shell, "intermediate/" + std::to_string(rank) + ".dat");
		outer_d_1_shell.clear();

		// Wait for all process to commit facets
		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0)
		{
			// Perform custom multifile sort on the intermediate simplexes also remove duplicate entries
			MultiFile<MapBinaryFile, std::pair<std::vector<short>, short>> facets("intermediate");
			facets.compressMap("input/" + std::to_string(iter_counter + 1) + ".dat", iter_counter);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		std::filesystem::remove("intermediate/" + std::to_string(rank) + ".dat");
		iter_counter++;
	}
#ifdef WRITE_CSV_OUTPUTS
	if (rank == 0)
	{
		MultiFile<VectorBinaryFile, std::vector<short>> dsimplexes("output");
		std::clog << "Found " << dsimplexes.writeCSV("dsimplexes.csv") << " dsimplexes for datset" << std::endl;
	}
#endif

	return;
}

// configPipe -> configure the function settings of this pipeline segment
template <typename nodeType>
bool helixDistPipe<nodeType>::configPipe(std::map<std::string, std::string> &configMap)
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
	this->ut.writeDebug("helixDistPipe", "Configured with parameters { eps: " + configMap["epsilon"] + " , debug: " + strDebug + ", outputFile: " + this->outputFile + " }");

	return true;
}
// outputData -> used for tracking each stage of the pipeline's data output without runtime
template <typename nodeType>
void helixDistPipe<nodeType>::outputData(pipePacket<nodeType> &inData)
{
	std::ofstream file;
	file.open("output/" + this->pipeType + "_output.csv");
	// code to print the data

	file.close();
	return;
}

template class helixDistPipe<simplexNode>;
template class helixDistPipe<alphaNode>;
template class helixDistPipe<witnessNode>;