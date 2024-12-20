/*
 * writeOutput hpp + cpp protoype and define a class for writing output from
 * the LHF system (https://github.com/wilseypa/LHF). The class is designed
 * to hold all functions corresponding to writing output data to .csv files,
 * .mat files, console, and other additional files as they are added to the system.
 *
 *
 *
 */

#include <string>
#include <iostream>
#include <vector>
#include <ctime>
#include <regex>
#include <fstream>
#include "writeOutput.hpp"

/**
 * @brief Construct a new write Output::write Output object
 *
 */
// writeOutput constructor, currently no needed information for the class constructor
writeOutput::writeOutput()
{
}

/**
 * @brief
 *
 * @param stats
 * @param filename
 * @return true
 * @return false
 */
// writeStat -> write the pipeline statistics to a csv formatted file
bool writeOutput::writeStats(const std::string &stats, const std::string &filename)
{
	std::ofstream file(filename + "_stats.csv");
	file << "PipeName,PhysExecutionTime,TransientMemorySize,MemoryUnits,VertexCount,SimplexCount\n";
	file << stats;
	file.close();
	return true;
}

/**
 * @brief
 *
 * @param stats
 * @param filename
 * @return true
 * @return false
 */
// writeRunLog -> write the run log to a csv formatted file
bool writeOutput::writeRunLog(const std::string &stats, const std::string &filename)
{
	std::ofstream file(filename + "_runLog.csv");
	file << "Timestamp,Process_#,Workload_#,PointSize,simplexCount,runtime(s),Arguments...\n";
	file << stats;
	file.close();
	return true;
}

/**
 * @brief
 *
 * @param data
 * @param filename
 * @return true
 * @return false
 */
bool writeOutput::writeBarcodes(const std::vector<bettiBoundaryTableEntry> &data, const std::string &filename)
{
	std::ofstream file(filename + ".csv");
	file << "dimension,birth,death\n";

	for (auto &row : data)
		file << std::to_string(row.bettiDim) << "," << std::to_string(row.birth) << "," << std::to_string(row.death) << '\n';

	file.close();
	return true;
}

/**
 * @brief
 *
 * @param data
 * @param filename
 * @return true
 * @return false
 */
// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(const std::vector<std::vector<double>> &data, const std::string &filename)
{
	std::ofstream file(filename, std::ios::app);
	if (!file.is_open())
	{
		std::cerr << "Failed to open the file for writing: " << filename << std::endl;
		return false;
	}

	for (auto &row : data)
	{
		for (size_t i = 0; i < row.size() - 1; ++i)
		{
			file << row[i] << ',';
		}
		if (!row.empty())
			file << row[row.size() - 1] << '\n';
	}
	file.close();
	return true;
}

/**
 * @brief
 *
 * @param data
 * @param filename
 * @return true
 * @return false
 */
// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(const std::string &data, const std::string &filename)
{
	std::ofstream file(filename + ".csv");
	file << data;
	file.close();
	return true;
}

/**
 * @brief
 *
 * @param data
 * @param filename
 * @param header
 * @return true
 * @return false
 */
// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(const std::vector<std::vector<double>> &data, const std::string &filename, const std::string &header)
{
	std::ofstream file(filename + ".csv");

	file << header;

	for (auto &row : data)
	{
		for (size_t i = 0; i < row.size() - 1; ++i)
		{
			file << row[i] << ',';
		}
		if (!row.empty())
			file << row[row.size() - 1] << '\n';
	}
	file.close();
	return true;
}

/**
 * @brief
 *
 * @param data
 * @param filename
 * @param header
 * @return true
 * @return false
 */
// writeCSV -> write a csv formatted file of data input
bool writeOutput::writeCSV(const std::string &data, const std::string &filename, const std::string &header)
{
	std::ofstream file(filename + ".csv");
	file << header;
	file << data;
	file.close();
	return true;
}

/**
 * @brief
 *
 * @param data
 * @param filename
 */
void writeOutput::writeBinary(const std::map<std::vector<short>, short> &data, const std::string &filename)
{
	if (data.empty())
		return;
	std::ofstream file(filename, std::ios::out | std::ios::binary);

	if (!file.is_open())
	{
		std::cerr << "Failed to open the file for writing." << std::endl;
		return;
	}

	// Write the size of the map
	const size_t mapSize = data.size();
	const size_t vectorSize = data.begin()->first.size();

	// Write the size of the map and vector
	file.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
	file.write(reinterpret_cast<const char *>(&vectorSize), sizeof(vectorSize));

	// Iterate through the map and write each key-value pair
	for (const auto &entry : data)
	{
		// Write the vector elements
		file.write(reinterpret_cast<const char *>(entry.first.data()), vectorSize * sizeof(short));
		// Write the short value associated with the vector
		file.write(reinterpret_cast<const char *>(&entry.second), sizeof(short));
	}

	file.close();
}

/**
 * @brief 
 * 
 * @param data 
 * @param filename 
 */
void writeOutput::writeBinary(const std::vector<std::pair<std::vector<short>, short>> &data, const std::string &filename)
{
	if (data.empty())
		return;
	std::ofstream file(filename, std::ios::out | std::ios::binary);

	if (!file.is_open())
	{
		std::cerr << "Failed to open the file for writing." << std::endl;
		return;
	}

	// Write the size of the map
	const size_t mapSize = data.size();
	const size_t vectorSize = data.begin()->first.size();

	// Write the size of the map and vector
	file.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
	file.write(reinterpret_cast<const char *>(&vectorSize), sizeof(vectorSize));

	// Iterate through the map and write each key-value pair
	for (const auto &entry : data)
	{
		// Write the vector elements
		file.write(reinterpret_cast<const char *>(entry.first.data()), vectorSize * sizeof(short));
		// Write the short value associated with the vector
		file.write(reinterpret_cast<const char *>(&entry.second), sizeof(short));
	}

	file.close();
}

/**
 * @brief
 *
 * @param data
 * @param filename
 */
void writeOutput::writeBinary(const std::vector<std::vector<short>> &data, const std::string &filename)
{
	if (data.empty())
		return;
	std::ofstream file(filename, std::ios::out | std::ios::binary);

	if (!file.is_open())
	{
		std::cerr << "Failed to open the file for writing." << std::endl;
		return;
	}

	// Write the size of the map
	const size_t dim1 = data.size();
	const size_t dim2 = data.begin()->size();

	// Write the size of the map and vector
	file.write(reinterpret_cast<const char *>(&dim1), sizeof(dim1));
	file.write(reinterpret_cast<const char *>(&dim2), sizeof(dim2));

	// Iterate through the map and write each key-value pair
	for (const auto &entry : data)
		file.write(reinterpret_cast<const char *>(entry.data()), dim2 * sizeof(short));

	file.close();
}

/**
 * @brief
 *
 * @param data
 * @param filename
 * @return true
 * @return false
 */
// writeMAT -> write in a mat formatted file of data input
bool writeOutput::writeMAT(const std::vector<std::vector<double>> &data, const std::string &filename)
{
	std::ofstream file(filename + ".mat");

	file << std::to_string(data.size()) << "\n"
		 << std::to_string(data[0].size()) << "\n";
	for (auto &row : data)
	{
		for (auto &column : row)
		{
			file << std::to_string(column) << "\n";
		}
	}
	file.close();
	return true;
}

/**
 * @brief
 *
 * @param data
 * @return true
 * @return false
 */
// writeConsole -> write data input to console (hopefully pretty-print)
bool writeOutput::writeConsole(const std::vector<bettiBoundaryTableEntry> &data)
{
	std::cout << "dimension,birth,death" << std::endl;

	for (auto &row : data)
		std::cout << std::to_string(row.bettiDim) << "," << std::to_string(row.birth) << "," << std::to_string(row.death) << std::endl;
	std::cout << std::endl;
	return true;
}

/**
 * @brief
 *
 * @param args
 * @param ident
 * @param wdStats
 * @param runtime
 * @return std::string
 */
std::string writeOutput::logRun(const std::map<std::string, std::string> &args, const std::string &ident, const std::string &wdStats, const std::string &runtime)
{
	// Get current time
	auto res = std::time(nullptr);

	std::string ret = std::to_string(res) + "," + ident + "," + wdStats + "," + runtime + ",";

	for (auto arg : args)
		ret += arg.first + ":" + arg.second + ",";

	ret += "\n";

	return ret;
}
