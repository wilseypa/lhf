# Python tools and scripts for testing and evaluating LHF

---

## CONTENTS

- DataScripts | Scripts for generating, formatting, and scaling data
  
- RunExperiments | Scripts for running LHF (and other tools) against data sets

- RunGraphs | Scripts for generating PDF result graphs from persistence intervals

- RunStats | Scripts for analysis of distance metrics between persistence intervals 

---
			  
### DataScripts 

- scaleColumns.py | Used to generate scaled and z-score normalized data from CSV source data

- tadasetsgen.py | Used to generate several TDA shapes defined in tadasets package

---

### RunExperiments

- runExperiments.py | Used to run different TDA libraries taking parameter inputs from shell script

- runMultiBatch.sh | Used to call runExperiments.py in a separate process; prevents test stoppage when python fails

---

### RunGraphs

- phOutputPlotting.py | Plots persistence intervals into persistence diagrams and barcode diagrams

---
 
### RunStats

- runUpscaleStats.py | Simplified python script for computing bottleneck, heat kernel, and wasserstein distances between persistence intervals

- runStats.sh | Used to call runAnalysis.py for computing distance metrics between persistence intervals

- runAnalysis.py | Probably unnecessary but runs analysis over multiple folders/files

---
