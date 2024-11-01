#pragma once

#include <vector>
#include <fstream>
#include <string>

struct VectorFile
{
    std::ifstream fileStream;
    std::vector<short> cache;

    VectorFile(const std::string &filename);
    bool operator<(const VectorFile &other) const;
    bool updateCache();
};

struct VectorBinaryFile
{
    std::ifstream fileStream;
    std::vector<short> cache;
    size_t dim1;
    size_t dim2;

    VectorBinaryFile(const std::string &filename);
    bool operator<(const VectorBinaryFile &other) const;
    bool updateCache();
};

struct MapBinaryFile
{
    std::pair<std::vector<short>, short> cache;
    std::ifstream fileStream;
    size_t mapSize;
    size_t vectorSize;

    MapBinaryFile(const std::string &filename);
    bool operator<(const MapBinaryFile &other) const;
    bool updateCache();
};

template <typename FileType, class baseType>
class MultiFile
{
private:
    std::vector<FileType> fileDataMap;
    baseType curr_element;
    bool readValue();
    bool readUnique();

public:
    MultiFile(const std::string &directory);
    void compressMap(const std::string &outputFileName, int iterationCounter);
    size_t writeCSV(const std::string &outputFileName);
};