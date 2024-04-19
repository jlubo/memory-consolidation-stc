# Memory recall dynamics for simulator comparison
***
####

The notebook __figures.ipynb__ simulates:
* learning of patterns of varied numbers of neurons with recall after 10s,
* learning of patterns of varied numbers of neurons, followed by consolidation and recall after 8h.

The model parameters and stimulation protocols are equal to those in Luboeinski & Tetzlaff, Commun. Biol., 2021 (https://doi.org/10.1038/s42003-021-01778-y). For each trial, a new network topology is generated.

The notebook here uses the C++ code and the build scripts found in the __simulation-code/__ folder of this repository, as well as the run scripts in the __simulation-bin/__ folder, and the script __averageFileColumnsAdvanced.py__ found in the __analysis/__ folder of this repository.

The results will be ready to be compared to results from comparable implementations, in particular, from the [Arbor simulator](https://github.com/jlubo/arbor_network_consolidation) or the [Brian 2 simulator](https://github.com/jlubo/brian_network_plasticity).

The comparison can be done using the scripts provided [here](https://github.com/jlubo/simulator_comparison).

