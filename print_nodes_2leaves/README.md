# Description:

The `print_nodes_2leaves.py` script reads in a user-specified unrooted
phylogenetic tree in the Newick format, midpoint roots the tree and searches
for nodes with two leaves. For each such node, prints out its leaves and the
corresponding bootstrap support.

# Requirements:

This script requires Python 3 (tested with Python 3.6.8 under Scientific Linux
release 6.8, Carbon) and [ETE toolkit](http://etetoolkit.org/).  This script
uses subroutine `search_by_size` from the ETE toolkit
[tutorial](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html).

# Usage:

`python print_nodes_2leaves.py input_nwk`

where:
`input_nwk` is the name of the input Newick formatted tree file


# Output files:

The script generates a single output text file with the extension
`.TWO.LEAVES.txt`.  Each line of the output file corresponds to a single node
with two leaves and contains names of these leaves and the corresponding
bootstrap support value.

# Sample input data:

A sample Newick formatted tree file (`sample.nwk`) is included along with the
script.

The command to run the script with the sample Newick file: 
`python print_nodes_2leaves.py sample.nwk`

The corresponding sample output file is provided in the subdirectory
`sample_output`.
