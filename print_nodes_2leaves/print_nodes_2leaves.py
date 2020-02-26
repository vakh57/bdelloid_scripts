#!/usr/bin/env python3

# This file is part of the software that was used to process sequencing data 
# for A. vaga individuals L1-L11 analyzed in the manuscript 
# https://www.biorxiv.org/content/10.1101/489393v1 
#
# Copyright (c) 2018-2020, by Olga Vakhrusheva <vakh57@gmail.com>
#
# The script is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# The script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with the software; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.

import sys

from ete3 import Tree

if len (sys.argv) != 2 :
    print ("\nError: no arguments provided or unexpected number of arguments (expecting a path to a file in the Newick format). Stopped.")
    sys.exit (1)

program_name = sys.argv[0]
arguments = sys.argv[1:]
nargs = len(arguments)

# print("nargs:", nargs)

NWK_TREE = sys.argv[1]

print("\n")
print("program_name:", program_name)

print("\n")
print("#############################################################################################")
print("Reads in a tree in the Newick format")
print("Midpoint roots the tree and looks for nodes with two leaves")
print("For each such node, prints out its leaves and the corresponding bootstrap support value")
print("#############################################################################################")

print("\n")
print("input Newick tree file:", NWK_TREE)

OUT_2LEAVES=NWK_TREE + '.TWO.LEAVES.txt'

print("\n")
print("output file:", OUT_2LEAVES)

OUT_2LEAVESF = open(OUT_2LEAVES,"w+")

### Reads a tree from the input Newick file
t = Tree(NWK_TREE)

### Midpoint roots the tree
MIDPR=t.get_midpoint_outgroup()

t.set_outgroup(MIDPR)


### This subroutine is from http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html
def search_by_size(node, size):
    "Finds nodes with a given number of leaves"
    matches = []
    for n in node.traverse("preorder"):
       if len(n) == size:
          matches.append(n)
    return matches


### Looks for nodes with 2 leaves
two_leaves=search_by_size(t, size=2)

### Iterates over all found nodes with 2 leaves, prints out the corresponding leaves and bootstrap support values
for node2l in two_leaves:
  # print("\n")
  # print("New node with 2 leaves:")
  # print(node2l.support)
  # print(node2l.get_leaves())

  node2lLeaves = []
  node2lLeaves=node2l.get_leaves()
  node2lLeaves_IDS = []

  for leaf in node2lLeaves:
    # print("leaf.name:",leaf.name)
    node2lLeaves_IDS.append(leaf.name)
  
  node2lLeaves_IDS_string = ' '.join(node2lLeaves_IDS)
  node2lLeaves_IDS_string = node2lLeaves_IDS_string + ";\t" + str(node2l.support)
  # print(node2lLeaves_IDS,node2l.support)   
  # print(node2lLeaves_IDS_string)   
  OUT_2LEAVESF.write(node2lLeaves_IDS_string + "\n") 

OUT_2LEAVESF.close() 

print("\n")
print("#############################################################################################")
print("Finished executing")
