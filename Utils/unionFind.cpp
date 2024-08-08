/**
 * @file unionFind.hpp
 *
 * @brief Contains unionFind class and supporting methods for the LHF system (https://github.com/wilseypa/LHF).
 */


#include <unionFind.hpp>

/**
 * @brief Initializes a Union-Find data structure.
 * 
 * @param n The number of elements in the Union-Find data structure.
 */
unionFind::unionFind(int n) : rank(n, 0), parent(n, 0)
{
	for (int i = 0; i < n; i++)
		parent[i] = i;
}

/**
 * @brief Finds the root of the component that contains the given element.
 * 
 * @param i The index of the element.
 * @return The index of the root of the component.
 */
int unionFind::find(int i)
{
	if (i == parent[i])
		return i;				 // Found name of the component
	parent[i] = find(parent[i]); // Path Compression
	return parent[i];
}

/**
 * @brief Joins two components together.
 * 
 * @param x The index of an element in one component.
 * @param y The index of an element in another component.
 * @return True if the two components were joined, false otherwise.
 * @return true 
 * @return false 
 */
bool unionFind::join(int x, int y)
{ // Union by rank
	x = find(x);
	y = find(y);
	if (x == y)
		return false;
	if (rank[x] == rank[y])
	{
		rank[y]++;
		parent[x] = y;
	}
	else if (rank[x] < rank[y])
	{
		parent[x] = y;
	}
	else
	{
		parent[y] = x;
	}
	return true;
}