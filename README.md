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
 |  "--dimensions" | "-d" | Max dimensions to run at | 3 |
 |  "--iterations" | "-r" | Number of iterations to run | 1000 |
 |  "--pipeline" | "-p" | Pipeline structure | default |
 |  "--inputFile" | "-i" | File to read into pipeline | None |
 |  "--outputFile" | "-o" | File to output to | None |

---
 
### EXAMPLES:

	##    ./LHF -p distMatrix --inputFile testData.csv
	##    ./LHF --pipeline distMatrix.distMatrix.distMatrix -i testData.csv -o output.csv

---
