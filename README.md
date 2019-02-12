# LHF: Lightweight Framework for Homology

---

### REQUIREMENTS 

- C++14
  
- CMake

---
			  
### COMPILING 

	##    cmake .
	##    make

---

###  RUNNING 

	##    cd LHF
	##    ./LHF <args>

---

### ARGUMENTS

 | Argument  | Shorthand | Description | Default
 | ------------- | ------------- | ------------- | ------------- |
 |  "--preprocessor" | "-pre" | Preprocessing Method | None |
 |  "--clusters" | "-k" | Number of preprocessing clusters | 5 |
 |  "--dimensions" | "-d" | Max dimensions to run at | 3 |
 |  "--epsilon" | "-e" | Epsilon value for simplicial complexes| 5 |
 |  "--complexType" | "-c" | Simplicial complex constructed| SimplexArrayList |
 |  "--upscaling" | "-u" | Upscaling selection (T/F) | True |
 |  "--iterations" | "-r" | Number of iterations to run | 1000 |
 |  "--pipeline" | "-p" | Pipeline structure | default |
 |  "--inputFile" | "-i" | File to read into pipeline | None |
 |  "--outputFile" | "-o" | File to output to | None |
 |  "--debug" | "-x" | Debug mode|0|

---
 
### EXAMPLES:

	##    ./LHF -p distMatrix --inputFile testData.csv
	##    ./LHF --pipeline distMatrix.distMatrix.distMatrix -i testData.csv -o output.csv

---
