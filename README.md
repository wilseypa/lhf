# LHF: Lightweight Framework for Homology

---

### REQUIREMENTS 

- C++17
  
- AutoTools

---
			  
### COMPILING 

	##    autoreconf -i
	##    ./configure
	##    make

---

###  RUNNING 

	##    cd src
	##    ./TDA_Cplusplus <args>

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

	##    ./TDA_Cplusplus -p distMatrix --inputFile testData.csv
	##    ./TDA_Cplusplus --pipeline distMatrix.distMatrix.distMatrix -i testData.csv -o output.csv

---
