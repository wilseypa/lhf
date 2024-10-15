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

```console
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
---

### PYTHON Package Local Build

```console
cp build/LHFmain/libLHFlib.so.1.0.0 pyLHF/src/lhf/libLHFlib.so
cd pyLHF/src/lhf
pip install .
```
---

###  RUNNING 

```console
cd LHFmain
./LHF <args>
```
---

### ARGUMENTS

 | Argument  | Shorthand | Default | Description | DataType
 | ------------- | ------------- | ------------- | ------------- | ------------- |
 | --betaMesh |-bmesh| null.csv| | |
 | --betaMode |-bm| noMode| | |
 | --beta |-b| 1| | `<int>`|
 | --alphaFilterationValue |-afv| 50000| | |
 | --nodeType |-n| simplexNode| | |
 | --reductionPercentage |-rp| 10| | |
 | --maxSize |-ms| 2000| | |
 | --threads |-t| 30| | `<int>`|
 | --threshold |-th| 250| | |
 | --scalar |-s| 0.5| | |
 | --mpi |-a| 0| | `<int>`|
 | --mode |-m| standard| Sets the mode for LHF to run in| (standard\|reduced\|upscale\|sw)|
 | --dimensions |-d| 1| Sets the maximum homology dimension to compute (H_d)| `<int>`|
 | --iterations |-r| 250| | `<int>`|
 | --pipeline |-p| | | |
 | --inputFile |-i| None| Filename (csv) for LHF input| `<filename>`|
 | --outputFile |-o| output| Filename for LHF output| `<filename>`|
 | --epsilon |-e| 5| Maximum epsilon threshold| `<float>`|
 | --lambda |-l| .25| Decay factor lambda for DenStream| |
 | --debug |-x| 0| | `<int(0\|1)>`|
 | --complexType |-c| simplexArrayList| | |
 | --clusters |-k| 20| | `<int>`|
 | --preprocessor |-pre| | | |
 | --upscale |-u| false| | `<bool>`|
 | --seed |-q| -1| | |
 | --twist |-w| false| | `<bool>`|
 | --collapse |-z| false| | `<bool>`|
 | --involutedUpscale |-iu| false| | `<bool>`|
 | --involuted| -inv| false| | `<bool>`|

---
 
### EXAMPLES:
```console
./LHF -m fast --inputFile testData.csv
./LHF --pipeline distMatrix.distMatrix.distMatrix -i testData.csv -o output.csv
```
---
