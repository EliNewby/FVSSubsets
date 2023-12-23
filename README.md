# FVSSubsets
This code was developed to identify subsets of the feedback vertex set (FVS) of a Boolean network and simulate how driving these subsets affects the dynamics of the network [1,2]

## Installation
This code functions best in an anaconda environment with these versions:
- Python 3.8.10
- Ginsim 3.0.0b (https://github.com/GINsim/GINsim-python)

and the FVS python code was developed by Gang Yang and is also available here: https://github.com/yanggangthu/FVS_python

## Functions
### FindBestSubsets
The FindBestSubsets file contains the function `findBestFVSSubsets` to find the topological metric values and intersection metric percentile cutoffs for FVS subsets of a specific size on a network in graphml format. The function returns a pandas DataFrame of each subset's topological values, a pandas DataFrame of each FVS subset's intersection metric percetile cutoffs, and a list of every FVS subset.

### RunSimulations
The RunSimulations file contains the function `calculateToandAwayControl` to simulate the network and calculate the *To Control* and *Away Control* values for a list of subsets. Along with a list of subsets, this code requires a Boolean network in one of the file types bioLQM accepts as input (see here: http://www.colomoto.org/biolqm/doc/formats.html). This function returns a pandas DataFrame of the *To Control* and *Away Control* values for each subset, and two dictionaries of Boolean values indicating if the subset is partially informative (goodAway) and if the subset is fully informative (goodTo).

### Example
The T-LGL 2 node subsets are run in the `FVSSubsets.ipynb` jupyter notebook file.
Along with running both of these functions, the notebook also recreates the figures from our paper.
The TLGLRandom folder contains two excel files of the *To Control* and *Away Control* values for the random samples of all subsets (`TLGLRandomResults.xlsx`) and FVS subsets (`TLGLRandomFVSResults.xlsx`) that are used to generate the box-and-whisker plots.

## Models
The models folder contains both the booleannet and graphml files for each of the models used.
The booleannet format is a Boolean network format which formats the Boolean rules as:
```
Node1 *= Node2 and Node3
Node2 *= Node1 or not Node3
Node3 *= not (Node1 and Node3)
```
where the text before the `*=` is the node name and the text after is the node's Boolean update rule using Boolean logical operations.

## References
[1] Eli Newby, Jorge Gómez Tejeda Zañudo, Réka Albert; Structure-based approach to identifying small sets of driver nodes in biological networks. *Chaos* 1 June 2022; 32 (6): 063102.

[2] Eli Newby, Jorge Gómez Tejeda Zañudo, Réka Albert; Structure-based approach to identify driver nodes in ensembles of biologically inspired Boolean networks. *Phys. Rev. Res.* Jul 2023; 5 (3):033009
