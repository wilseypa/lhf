#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <iostream>
#include <algorithm>
#include "multifileops.hpp"

/**
 * @brief Construct a new Vector File:: Vector File object
 *
 * @param filename
 */
VectorFile::VectorFile(const std::string &filename) : fileStream(filename, std::ios::in)
{
    if (!fileStream.is_open())
        throw std::runtime_error("Failed to open the file.");
    updateCache();
};

/**
 * @brief
 *
 * @param other
 * @return true
 * @return false
 */
bool VectorFile::operator<(const VectorFile &other) const
{
    return cache < other.cache;
};

/**
 * @brief
 *
 * @return true
 * @return false
 */
bool VectorFile::updateCache()
{
    std::string line;
    if (std::getline(fileStream, line))
    {
        cache.clear();
        std::string cell;
        std::istringstream lineStream(line);
        while (std::getline(lineStream, cell, ','))
            cache.push_back(std::stoi(cell));
        return true;
    }
    else
    {
        fileStream.close();
        return false;
    }
};

/**
 * @brief Construct a new Vector Binary File:: Vector Binary File object
 *
 * @param filename
 */
VectorBinaryFile::VectorBinaryFile(const std::string &filename) : fileStream(filename, std::ios::in | std::ios::binary)
{
    if (!fileStream.is_open())
        throw std::runtime_error("Failed to open the file.");
    fileStream.read(reinterpret_cast<char *>(&dim1), sizeof(dim1));
    fileStream.read(reinterpret_cast<char *>(&dim2), sizeof(dim2));
    cache.resize(dim2);
    updateCache();
};

/**
 * @brief
 *
 * @param other
 * @return true
 * @return false
 */
bool VectorBinaryFile::operator<(const VectorBinaryFile &other) const
{
    return cache < other.cache;
};

/**
 * @brief
 *
 * @return true
 * @return false
 */
bool VectorBinaryFile::updateCache()
{
    if (dim1 > 0)
    {
        fileStream.read(reinterpret_cast<char *>(cache.data()), dim2 * sizeof(short));
        dim1--;
        return true;
    }
    else
    {
        fileStream.close();
        return false;
    }
};

/**
 * @brief Construct a new Map Binary File:: Map Binary File object
 *
 * @param filename
 */
MapBinaryFile::MapBinaryFile(const std::string &filename) : fileStream(filename, std::ios::in | std::ios::binary)
{
    if (!fileStream.is_open())
        throw std::runtime_error("Failed to open the file." + filename);
    fileStream.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    fileStream.read(reinterpret_cast<char *>(&vectorSize), sizeof(vectorSize));
    cache.first.resize(vectorSize);
    updateCache();
};

/**
 * @brief
 *
 * @param other
 * @return true
 * @return false
 */
bool MapBinaryFile::operator<(const MapBinaryFile &other) const
{
    return cache.first < other.cache.first;
};

/**
 * @brief
 *
 * @return true
 * @return false
 */
bool MapBinaryFile::updateCache()
{
    if (mapSize > 0)
    {
        fileStream.read(reinterpret_cast<char *>(cache.first.data()), vectorSize * sizeof(short));
        fileStream.read(reinterpret_cast<char *>(&cache.second), sizeof(short));
        mapSize--;
        return true;
    }
    else
    {
        fileStream.close();
        return false;
    }
};

/**
 * @brief
 *
 * @tparam FileType
 * @tparam baseType
 * @return true
 * @return false
 */
template <typename FileType, class baseType>
bool MultiFile<FileType, baseType>::readValue()
{
    if (fileDataMap.empty())
        return false;
    auto minElementIt = std::min_element(fileDataMap.begin(), fileDataMap.end());
    curr_element = minElementIt->cache;
    if (!minElementIt->updateCache())
        fileDataMap.erase(minElementIt);
    return true;
};

/**
 * @brief
 *
 * @tparam FileType
 * @tparam baseType
 * @return true
 * @return false
 */
template <typename FileType, class baseType>
bool MultiFile<FileType, baseType>::readUnique()
{
    while (readValue())
    {
        bool isUnique = true;
        for (auto it = fileDataMap.begin(); it != fileDataMap.end(); ++it)
        {
            if (curr_element.first == it->cache.first)
            {
                isUnique &= curr_element.second == it->cache.second;
                if (!it->updateCache())
                    fileDataMap.erase(it--);
            }
        }
        if (isUnique)
            return true;
    }
    return false;
};

/**
 * @brief Construct a new Multi F Ile< File Type, base Type>:: Multi F Ile object
 *
 * @tparam FileType
 * @tparam baseType
 * @param directory
 */
template <typename FileType, class baseType>
MultiFile<FileType, baseType>::MultiFile(const std::string &directory)
{
    if (!std::filesystem::is_directory(directory))
        throw std::invalid_argument("Not a valid directory: " + directory);
    for (const auto &entry : std::filesystem::directory_iterator(directory))
        if (entry.is_regular_file())
            fileDataMap.emplace_back(entry.path().string());
};

/**
 * @brief
 *
 * @tparam FileType
 * @tparam baseType
 * @param outputFileName
 * @param iterationCounter
 */
template <typename FileType, class baseType>
void MultiFile<FileType, baseType>::compressMap(const std::string &outputFileName, int iterationCounter)
{
    // Open the output file for writing
    std::ofstream outputFile(outputFileName, std::ios::out | std::ios::binary);

    // Read from the previous iteration's data file
    FileType previousReader("input/" + std::to_string(iterationCounter) + ".dat");

    // Initialize variables for map size and vector size
    size_t mapSize = 0, vectorSize = previousReader.cache.first.size();

    outputFile.write(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));       // Write a placeholder for map size to the output file
    outputFile.write(reinterpret_cast<char *>(&vectorSize), sizeof(vectorSize)); // Write the size of the vector to facilitate reconstruction

    // Read values and write them to the output file
    while (readUnique())
    {
        if (curr_element.second != -1)
        {
            while (previousReader.cache.first < curr_element.first && previousReader.updateCache())
                continue;
            // Seek the reader from the previous iteration to the current value
            if (previousReader.cache.first != curr_element.first) // Write binary to the file if curr_element was not processed in the previous iteration
            {
                outputFile.write(reinterpret_cast<const char *>(curr_element.first.data()), vectorSize * sizeof(short));
                outputFile.write(reinterpret_cast<const char *>(&curr_element.second), sizeof(short));
                mapSize++;
            }
        }
    }

    // Seek to the beginning of the file and overwrite the mapSize binary with the original value
    outputFile.seekp(0);
    outputFile.write(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    outputFile.close();

    // Close the file streams and remove the previous iteration's data file
    previousReader.fileStream.close();
    std::filesystem::remove("input/" + std::to_string(iterationCounter) + ".dat");
};

/**
 * @brief
 *
 * @tparam FileType
 * @tparam baseType
 * @param outputFileName
 * @return size_t
 */
template <typename FileType, class baseType>
size_t MultiFile<FileType, baseType>::writeCSV(const std::string &outputFileName)
{
    std::ofstream outputFile(outputFileName, std::ios::out);
    size_t size = 0;
    while (readValue())
    {
        for (auto it = fileDataMap.begin(); it != fileDataMap.end(); ++it)
        {
            if (curr_element == it->cache && !it->updateCache())
                fileDataMap.erase(it--);
        }
        for (size_t i = 0; i < curr_element.size() - 1; ++i)
            outputFile << curr_element[i] << ',';                    // Add a comma to separate values (CSV format)
        outputFile << curr_element[curr_element.size() - 1] << '\n'; // Add a newline character to separate rows
        size++;
    }
    outputFile.close();
    return size;
};