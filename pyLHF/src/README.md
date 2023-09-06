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

`pip install lhf`

---
 
## Documentation:

LHF docs are included with the source repository and built with doxygen. To build and view the docs, clone the 
[LHF github repository](http//github.com/wilseypa/lhf). Build the documentation from the root folder using:

`doxygen`
    
Doxygen creates both HTML and Latex documentation in the *lhf/docs* folder. For HTML, navigate to index.html and open in a web browser
to view the class documentation of both the python interface and the LHF cpp library. 


---
 
## Examples:

### Run PH on input data

```
   import lhf
   import tadasets

   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
   pyLHF = lhf.createPipeline()

   #Set LHF args
   pyLHF.config["epsilon"] = 1.0
   pyLHF.config["dimensions"] = 3

   boundaryPis, pipePacket, elapsed = pyLHF.runPH(data)
```

### Plot the output in several ways

```
   import lhf
   import tadasets
   import numpy as np

   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
   pyLHF = lhf.createPipeline()

   boundaryPis, pipePacket, elapsed = pyLHF.runPH(data)

   #Remove boundary labels from the PIs
   pis = np.array([[z[0], z[1], z[2]] for z in boundaryPis])

   LHF.OutputAnalysis.persistenceDiagram(pis)
   LHF.OutputAnalysis.barcodeDiagrams(pis)
   LHF.OutputAnalysis.bettiCurve(pis)
```

### Run Partitioned Persistent Homology

```
   import lhf
   import tadasets

   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
   pyLHF = lhf.createPipeline()

   #Set LHF args
   pyLHF.config["epsilon"] = 1.0
   pyLHF.config["dimensions"] = 3
   pyLHF.config["mode"] = "iterUpscale"
   pyLHF.config["preprocessor"] = "kmeans++"
   pyLHF.config["clusters"] = 20

   boundaryPis, pipePacket, elapsed = pyLHF.runPH(data)
```

### Run PH on Weighted Alpha Complex

```
   import lhf
   import tadasets

   data = tadasets.dsphere(n=50, d=5, r=1, noise=0.1) 
   pyLHF = lhf.createPipeline()

   #Set LHF args
   pyLHF.config["epsilon"] = 1.0
   pyLHF.config["dimensions"] = 3
   pyLHF.config["mode"] = "weightedAlpha"

   boundaryPis, pipePacket, elapsed = pyLHF.runPH(data)
```

---
