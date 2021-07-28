# FVSSubsets
## Installation
This code functions best in an anaconda environment with these versions:
- Python 3.9.6
- Ginsim 3.0.0b (https://anaconda.org/colomoto/ginsim)

and the FVS python code is available here: https://github.com/yanggangthu/FVS_python

The FindBestSubsets file contains the code to find the topological metric values and intersection percentile cutoffs for a network. To run it, a network in graphml format is necessary.

The RunSimulations file contains the code to simulate the network and calculate the *To Control* and *Away Control* values for a list of subsets. Along with a list of subsets, this code requires a Boolean network input in one of the file types bioLQM accepts as input (see here: http://www.colomoto.org/biolqm/doc/formats.html)

For more information and an example, see the jupyter notebook file.
