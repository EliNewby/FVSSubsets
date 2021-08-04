# FVSSubsets
## Installation
This code functions best in an anaconda environment with these versions:
- Python 3.8.10
- Ginsim 3.0.0b (https://github.com/GINsim/GINsim-python)

and the FVS python code is available here: https://github.com/yanggangthu/FVS_python

### FindBestSubsets
The FindBestSubsets file contains the code to find the topological metric values and intersection percentile cutoffs for FVS subsets of a specific size on a network in graphml format. The function has an additional parameter of a list of intersections that are of interest, where each intersection is a list of metrics. The default parameters test the intersection of all seven metrics, the intersection of the propagation metrics, and the intersection of the Modified PRINCE and CheiRank metrics. To test other intersections, the metrics must be written as `OutDegree`, `Distance`, `PRINCE`, `ModPRINCE`, `CheiRank`, `Cycles`, and `SCC` to be recognized.

The function returns a pandas DataFrame of each subset's topological values, a pandas DataFrame of each subset's intersection percetile cutoffs, and a list of every FVS subset.

### RunSimulations
The RunSimulations file contains the code to simulate the network and calculate the *To Control* and *Away Control* values for a list of subsets. Along with a list of subsets, this code requires a Boolean network in one of the file types bioLQM accepts as input (see here: http://www.colomoto.org/biolqm/doc/formats.html). There are two more parameters in the function, which determine the amout of simulations to be run. The first parameter `wtSims` indicates the number of simulations to run on the wild-type system to get the basins of attraction of the network's attractors (1000 in our paper). The second parameter `subsetSims` indicates the number of simulations to run on each intervention (100 in our paper).

This function returns a pandas DataFrame of the *To Control* and *Away Control* values for each subset, and two dictionaries of Boolean values indicating if the subset doesn't drive to all of the attractors (goodAway) and if the subset drives to any (but not all) of the attractors (goodTo).

### Example
The T-LGL 2 node subsets are run in the `FVSSubsets.ipynb` jupyter notebook file. Along with running both of these functions, the notebook also recreates the figures from our paper.

## Models
The models folder contains both the booleannet and graphml files for each of the models used.
The booleannet format is a Boolean network format which formats the Boolean rules as:
```
Node1 *= Node2 and Node3
Node2 *= Node1 or not Node3
Node3 *= not (Node1 and Node3)
```
where the text before the `*=` is the node name and the text after is the node's Boolean update rule using Boolean logical operations.
