#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "witnessComplex.hpp"

// Yanxio , Yanlin and  Jisun are working on witness complex and will post there code here.

template <typename nodeType>
witnessComplex<nodeType>::witnessComplex(double maxE, double maxD) : simplexArrayList<nodeType>::simplexArrayList(0, 0)
{
	/**
	witnessComplex(double maxE, double maxD)

	@brief Initializes the witness complex
	@tparam nodeType The data type of the simplex node.
	@param maxE The max epsilon limit for complex construction.
	@param maxD The max dimension limit for complex construction.
*/

	std::cout << "Constructed Witness Complex!" << std::endl;

	this->simplexType = "witnessComplex";
	this->maxEpsilon = maxE;
	this->maxDimension = maxD;
}

template <typename nodeType>
witnessComplex<nodeType>::~witnessComplex()
{
	/**
	~witnessComplex()

	@brief Destructs the witness complex and frees used memory.
	@tparam nodeType The data type of the simplex node.
*/
	this->simplexList.clear();
}

// Explicit Template Class Instantiation
template class witnessComplex<simplexNode>;
template class witnessComplex<alphaNode>;
template class witnessComplex<witnessNode>;
