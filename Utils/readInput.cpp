/*
 * readInput hpp + cpp protoype and define a class for reading input into
 * the LHF system (https://github.com/wilseypa/LHF). The class is designed
 * to hold all functions corresponding to reading input data from .csv files,
 * .mat files, and other additional files as they are added to the system.
 *
 *
 *
 */

/**
  @file readInput.hpp
  @brief Defines a class for reading input into the LHF system (https://github.com/wilseypa/LHF).
 */

#include "readInput.hpp"

/**
  @brief Constructor for readInput class.

  Currently no needed information for the class constructor.
 */

readInput::readInput()
{
}

/**
  @brief Reads in a csv formatted file from file input.

  @param filename The complete filename and relative path (if needed) for reading.
  @return A vector array of doubles.

  The file format is:
  @code
  1.0423, 1.0244, 1.032 \n
  @endcode
  By default, double conversion (std::stod) will handle 'E' and 'e' for scientific notation.

  Strips whitespace characters where appropriate to create a vector array of doubles.

  @todo Handle if there is a comma at the end of the line (i.e. "5,5,5," should be <5,5,5>).
 */

std::vector<std::vector<double>> readInput::readCSV(const std::string &filename)
{
	std::vector<std::vector<double>> result;

	std::ifstream file(filename);

	if (!file)
	{
		std::cout << "Failed to open file: " << filename << std::endl;
		return result;
	}
	size_t vec_len = 0;
	// We are going to iterate through each line of the file until we reach the end
	while (!file.eof())
	{
		std::string line;		 // Temporary (current) line
		std::vector<double> tmp; // Temporary (current) vector
		getline(file, line);	 // Read the next line from file
		if (parseDoubleVector(line, tmp))
		{
			if (vec_len == tmp.size() || vec_len == 0, vec_len = tmp.size())
				result.push_back(tmp);
			else
			{
				std::cout << "Vectors in CSV are expected to be of same length\n"
						  << "First vector was of size " << vec_len << ", but found vector of size " << tmp.size() << " in file." << std::endl;
				break;
			}
		}
	}
	file.close();

	return result;
}

/**
  @brief Parses a double vector from a given line of text.

  @param line The line of text to parse.
  @param row A reference to the output vector array of doubles.
  @return True if the vector was parsed successfully, false otherwise.
 */

bool readInput::parseDoubleVector(const std::string & readline, std::vector<double> &row)
{
	std::size_t pos = std::string::npos;

	// Replace whitespace in the current line
	std::string line = std::regex_replace(readline, std::regex(" "), "");

	// Check if the line has a length (is not a blank line)
	if (line.size() > 1)
	{

		// Iterate through each comma of the csv
		while ((pos = line.find_first_of(",")) != std::string::npos)
		{
			// Push the value found before the comma, remove from the line
			row.push_back(std::stod(line.substr(0, pos)));
			line.erase(0, pos + 1);
		}
		// Get the last value of the line
		if (line.size() > 0)
			row.push_back(std::stod(line.substr(0, pos)));
	}
	else
		return false;

	return true;
}

/**
  @brief Reads in a mat formatted file from file input.

  @param filename The complete filename and relative path (if needed) for reading.
  @return A vector array of doubles.

  The file format is:
  @code
  15		(# of vectors)
  20		(# of dimensions)
  1.023	(value of [0,0])
  4.234	(value of [0,1])
  ...	(subsequent values of [i, j])
  @endcode
  By default, double conversion (std::stod) will handle 'E' and 'e' for scientific notation.

  Strips whitespace characters where appropriate to create a vector array of doubles.

  @todo A lot.
 */

std::vector<std::vector<double>> readInput::readMAT(const std::string &filename)
{
	std::vector<std::vector<double>> result;
	int vectors = 0;
	int dimensions = 0;

	std::ifstream file;
	file.open(filename);

	std::string line; // Temporary (current) line

	// Get the number of vectors
	if (getline(file, line))
	{
		line = std::regex_replace(line, std::regex(" "), "");
		vectors = std::stoi(line);
	}
	else
	{
		std::cout << "Expected file to begin with vector length" << std::endl;
		file.close();
		return result;
	}
	// Get the number of dimensions
	if (getline(file, line))
	{
		line = std::regex_replace(line, std::regex(" "), "");
		dimensions = std::stoi(line);
	}
	else
	{
		std::cout << "Dimensions of matrix expected to present in file" << std::endl;
		file.close();
		return result;
	}

	// We are going to iterate through each line of the file until we reach the end
	for (int vect = 0; vect < vectors; vect++)
	{
		std::vector<double> tmp; // Temporary (current) vector

		for (int dim = 0; dim < dimensions; dim++)
		{
			getline(file, line); // Read the next line from file

			// Replace whitespace in the current line
			line = std::regex_replace(line, std::regex(" "), "");

			// Check if the line has a length (is not a blank line)
			if (line.size() > 0)
			{

				// Push the value found before the comma, remove from the line
				tmp.push_back(std::stod(line));
			}
		}
		result.push_back(tmp);
	}
	file.close();
	return result;
}

/**
  @brief Initializes the input stream for reading from a file.

  @param filename The complete filename and relative path (if needed) for reading.
  @return True if the stream was initialized successfully, false otherwise.
 */

bool readInput::streamInit(const std::string &filename)
{
	pFile = fopen(filename.c_str(), "r");

	if (pFile == NULL)
		return false;

	return true;
}

/**
  @brief Reads the next row of data from the input stream.

  @param row A reference to the output vector array of doubles.
  @return True if the row was read successfully, false otherwise.
 */

bool readInput::streamRead(std::vector<double> &row)
{
	if (pFile == NULL)
		return false;

	// Check for EOF
	if (fgets(streamBuffer, 1000, pFile) == NULL)
		return false;

	// Parse our row
	if (parseDoubleVector(streamBuffer, row))
		return true;
	else
		return false;
}
