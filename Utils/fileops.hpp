#include <string>
#include <map>
#include <set>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iostream>

void writeCSV(const std::vector<std::vector<short>> &data, const std::string &filename)
{
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open())
    {
        std::cerr << "Failed to open the file for writing: " << filename << std::endl;
        return;
    }

    for (const std::vector<short> &row : data)
    {
        for (size_t i = 0; i < row.size() - 1; ++i)
        {
            file << row[i] << ','; // Add a comma to separate values (CSV format)
        }
        file << row[row.size() - 1] << '\n'; // Add a newline character to separate rows
    }
    file.close();
}

std::vector<std::vector<short>> readVectorCSV(const std::string &fileName)
{
    std::vector<std::vector<short>> data;       // This will store the CSV data
    std::ifstream file(fileName, std::ios::in); // Open the CSV file
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return data; // Return an empty map if the file cannot be opened
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::vector<short> key;
        std::istringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ','))
            key.push_back((short)std::stoi(cell)); // Add to the key vector
        data.emplace_back(key);
    }
    file.close(); // Close the CSV file
    return data;
}

void writeBinaryFile(const std::map<std::vector<short>, short> &data, const std::string &filename)
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

void writeBinaryFile(const std::vector<std::vector<short>> &data, const std::string &filename)
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

std::map<std::vector<short>, short> readBinaryMapFile(const std::string &filename, int numProcesses = 1, int rank = 0)
{
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Failed to open the file for reading." << std::endl;
        return {};
    }

    std::map<std::vector<short>, short> data;

    size_t mapSize;
    size_t vectorSize;

    // Read the size of the map and vector
    file.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    file.read(reinterpret_cast<char *>(&vectorSize), sizeof(vectorSize));
    size_t block_size = ceil((double) mapSize / numProcesses);
    size_t end = block_size * (rank+1);
    size_t start = (rank == 0) ? 0 : end - block_size;

    // Read data using standard file operations
    std::vector<short> key(vectorSize);
    short value;

    // Seek to the correct position in the file
    file.seekg(sizeof(size_t) * 2 + (start * (vectorSize + 1)) * sizeof(short));

    for (size_t i = start; i < end; ++i)
    {
        file.read(reinterpret_cast<char *>(key.data()), vectorSize * sizeof(short));
        file.read(reinterpret_cast<char *>(&value), sizeof(short));
        data.emplace(key, value);
    }

    file.close();
    return data;
}

std::vector<std::vector<short>> readBinaryVectorFile(const std::string &filename)
{
    std::ifstream file(filename, std::ios::in | std::ios::binary);

    if (!file.is_open())
    {
        std::cerr << "Failed to open the file for reading." << std::endl;
        return {};
    }

    size_t dim1;
    size_t dim2;

    // Read the dimensions of the vector
    file.read(reinterpret_cast<char *>(&dim1), sizeof(dim1));
    file.read(reinterpret_cast<char *>(&dim2), sizeof(dim2));

    std::vector<std::vector<short>> data;
    data.reserve(dim1);

    while (dim1-- > 0)
    {
        std::vector<short> vec(dim2);
        file.read(reinterpret_cast<char *>(vec.data()), dim2 * sizeof(short));
        data.emplace_back(vec);
    }

    file.close();
    return data;
}