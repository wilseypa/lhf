
/* 
 * Need for data structure for d-way space partioning:
 * Issues with current partioning trees:
 * 1. Existing data structures have worst time complexity in higher dimensions.
 * 2. Data points required explods with dimension 2^d.
 * 3. Have same time complexity as bruit force in higher dimensions.
 * 4. Partioning of space increases with dimensions 2^d
 * 
 * How d-way trees can be fruitful in higher dimensionsl space.
 * 1. The space partioning is leaner to the data diension.
 * 2. Required point is O(d+1) as compared to O(2^d).
 * 3. Is a way to provide standary multi-way partioning.
 * 4. NN search complexity can be decreased.
 */ 


/*
				*************************Polytopal Complex***************************
*1. Will be useful in convex partioning of space.
*2. Could be an alternative/assist in faster approximate Delaunay Triangulation.
*3. Partion space in convex components.
*4. Highly Scalable with NN searches.

*/
