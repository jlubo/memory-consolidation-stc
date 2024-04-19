# Basic dynamics for simulator comparison
***
####

The notebook __figures.ipynb__ simulates:
* basic early-phase dynamics,
* basic late-phase plasticity dynamics,
* `smallnet3` dynamics.

The model parameters are equal to those in Luboeinski & Tetzlaff, Commun. Biol., 2021 (https://doi.org/10.1038/s42003-021-01778-y).

The notebook here uses the C++ code and the build scripts found in the __simulation-code/__ folder of this repository, as well as the run scripts in the __simulation-bin/__ folder, and the script __extractAverageQMI.py__ found in the __analysis/__ folder of this repository.

The results will be ready to be compared to results from comparable implementations, in particular, from the [Arbor simulator](https://github.com/jlubo/arbor_network_consolidation) or the [Brian 2 simulator](https://github.com/jlubo/brian_synaptic_plasticity_stc). 

The comparison can be done using the scripts provided [here](https://github.com/jlubo/simulator_comparison).

