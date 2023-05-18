# LHF: Lightweight Homology Framework

LHF is a homology framework designed for enabling modular pipelines for preprocessing and experiments around computing persistent homology. 
The base pipelines enable LHF to compute persistence intervals of an input point cloud by evaluating the distance matrix, building the 
simplicial complex, and reducing the complex to identify the persistence intervals at different dimensions. 

Additional pipelines have been created for preprocessing, approximations of the boundary matrix, approximate simplicial complexes 
representing the input point cloud, boundary extraction, and upscaling based on previous studies detailed in several of the references 
below. The intent for LHF is to provide a modular framework to continue exploring, developing, and evaluating mechanisms for approximating 
or optimizing the computation of persistent homology over an input point cloud.

The python interface to LHF, available from [pypi](http://pypi.org/project/lhf/), uses ctypes to expose the LHF shared object library. This enables python users
to call the LHF library on data and retrieve the entire pipepacket object used by LHF. 

---
  
## Installing:

$ pip install lhf

---
 
## Documentation:

LHF docs are included with the source repository and built with doxygen. To build and view the docs, clone the 
[LHF github repository](http//github.com/wilseypa/lhf). Build the documentation from the root folder using:

$   doxygen
    
Doxygen creates both HTML and Latex documentation in the *lhf/docs* folder. For HTML, navigate to index.html and open in a web browser
to view the class documentation of both the python interface and the LHF cpp library. 


---
 
## Examples:

### Run PH on input data

$   from LHF import LHF
$   import tadasets
$
$   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
$   pyLHF = LHF.LHF()
$
$   #Set LHF args
$   pyLHF.args["epsilon"] = 1.0
$   pyLHF.args["dimensions"] = 3
$
$   boundaryPis, pipePacket = pyLHF.runPH(data)

### Plot the output in several ways

$   from LHF import LHF
$   from LHF.OutputAnalysis import persistenceDiagram, barcodeDiagram, bettiCurve
$   import tadasets
$   import numpy as np
$
$   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
$   pyLHF = LHF.LHF()
$
$   boundaryPis, pipePacket = pyLHF.runPH(data)
$
$   #Remove boundary labels from the PIs
$   pis = np.array([[z[0], z[1], z[2]] for z in boundaryPis])
$
$   persistenceDiagram(pis)
$   barcodeDiagrams(pis)
$   bettiCurve(pis)

### Run Partitioned Persistent Homology

$   from LHF import LHF
$   import tadasets
$
$   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
$   pyLHF = LHF.LHF()
$
$   #Set LHF args
$   pyLHF.args["epsilon"] = 1.0
$   pyLHF.args["dimensions"] = 3
$   pyLHF.args["mode"] = "iterUpscale"
$   pyLHF.args["preprocessor"] = "kmeans++"
$   pyLHF.args["clusters"] = 20
$
$   boundaryPis, pipePacket = pyLHF.runPH(data)

### Run PH on Weighted Alpha Complex

$   from LHF import LHF
$   import tadasets
$
$   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
$   pyLHF = LHF.LHF()
$
$   #Set LHF args
$   pyLHF.args["epsilon"] = 1.0
$   pyLHF.args["dimensions"] = 3
$   pyLHF.args["mode"] = "weightedAlpha"
$
$   boundaryPis, pipePacket = pyLHF.runPH(data)

---
