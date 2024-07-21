# LHF: Lightweight Homology Framework

LHF is a homology framework designed for enabling modular pipelines for preprocessing and experiments around computing persistent homology. The base pipelines enable LHF to compute persistence intervals of an input point cloud by evaluating the distance matrix, building the simplicial complex, and reducing the complex to identify the persistence intervals at different dimensions. 

Additional pipelines have been created for preprocessing, approximations of the boundary matrix, approximate simplicial complexes representing the input point cloud, boundary extraction, and upscaling based on previous studies detailed in several of the references below. The intent for LHF is to provide a modular framework to continue exploring, developing, and evaluating mechanisms for approximating or optimizing the computation of persistent homology over an input point cloud.

---

### REQUIREMENTS 

- C++23
  
- CMake

- OpenMP

- Eigen3

- MPI (MPICH, OpenMPI)

- CGAL

---
			  
### COMPILING 

	mkdir build
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=Release
	make

---

###  RUNNING 

	cd LHFmain
	./LHF <args>

---

### ARGUMENTS

 | Argument  | Shorthand | Description | Default
 | ------------- | ------------- | ------------- | ------------- |
 |  "--preprocessor" | "-pre" | Preprocessing Method | None |
 |  "--clusters" | "-k" | Number of preprocessing clusters | 5 |
 |  "--dimensions" | "-d" | Max dimensions to run at | 3 |
 |  "--epsilon" | "-e" | Epsilon value for simplicial complexes| 5 |
 |  "--lambda" | "-l" | Lambda value (decay factor) for DenStream | .25 |
 |  "--mode" | "-m" | Mode to run LHF in (fast, slidingwindow, upscaling, etc.) | default |
 |  "--complexType" | "-c" | Simplicial complex constructed| SimplexArrayList |
 |  "--inputFile" | "-i" | File to read into pipeline | None |
 |  "--outputFile" | "-o" | File to output to | None |
 |  "--debug" | "-x" | Debug mode|0|

---
 
### EXAMPLES:

	##    ./LHF -m fast --inputFile testData.csv
	##    ./LHF --pipeline distMatrix.distMatrix.distMatrix -i testData.csv -o output.csv

---
