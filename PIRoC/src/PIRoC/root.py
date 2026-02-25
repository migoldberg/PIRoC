"""
root.py
Attempts to root the tree at an outgroup if possible.
"""

from inspect import ismodule
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree
from .metadata import parse_species_name, get_group_from_species_name

def find_monophyletic_outgroup(tree, outgroup_leaves):
    """
    Tries to find the most monophyletic clade of outgroup sequences in the tree.
    """

    # creates a set of the outgroup leaves
    outgroup_set = set(outgroup_leaves)

    # the most monophyletic outgroup node
    monophyletic_outgroup_node = None

    # the number of sequences in the most monophyletic outgroup node
    monophletic_outgroup_node_leaves = 0

    # loops throug each node in the tree in postorder
    for idx, node in enumerate(tree.traverse("postorder")):
        # adds an index feature to the node
        node.add_feature("index", idx)

        if node.is_leaf():
            # if the node is a leaf then there is one leaf
            num_leaves = 1

            # if the node is an outgroup leaf then there is one outgroup leaf
            num_outgroup_leaves = 1 if node in outgroup_set else 0
        else:
            # if the node is not a leaf then there are multiple leaves
            num_leaves = sum(child.num_leaves for child in node.children)

            # finds the number of outgroup leaves which are children of the node
            num_outgroup_leaves = sum(child.num_outgroup_leaves for child in node.children)

        if node.up is not None:
            # if the node has a parent then add the number of leaves and outgroup leaves to the node features
            node.add_feature("num_leaves", num_leaves)
            node.add_feature("num_outgroup_leaves", num_outgroup_leaves)

            # calculates the fraction of leaves in the node that are outgroups
            outgroup_fraction = (num_outgroup_leaves / num_leaves) if num_leaves else 0.0

            # adds the outgroup fraction to the node features
            node.add_feature("outgroup_fraction", outgroup_fraction)

            # checks if the node is monophyletic
            is_monophyletic = (num_outgroup_leaves == num_leaves) and (num_outgroup_leaves > 0)

            # if the node is more monophyletic than the current most monophyletic outgroup node then update the most monophyletic outgroup node
            if is_monophyletic and num_outgroup_leaves > monophletic_outgroup_node_leaves:
                monophletic_outgroup_node_leaves = num_outgroup_leaves
                monophyletic_outgroup_node = node

    # returns the most monophyletic outgroup node
    return monophyletic_outgroup_node

def root_tree_at_outgroup(tree, outgroup_leaves):
    """
    Attempts to root the tree at a monophyletic outgroup.
    """

    # assumes the tree is not rooted (by this script, if this does not re-root it than the tree-builder may have rooted it somewhere)
    rooted = False

    # if there are no outgroup leaves then return the tree
    if not outgroup_leaves:
        return tree, rooted

    # if there is only one outgroup leaf then root the tree there and return the tree
    if len(outgroup_leaves) == 1:
        outgroup_leaf = outgroup_leaves[0]
        tree.set_outgroup(outgroup_leaf)
        rooted = True
        return tree, rooted

    try:
        # finds the most monophyletic outgroup
        monopheltic_outgroup = find_monophyletic_outgroup(tree, outgroup_leaves)

        # root the tree at the most monophyletic outgroup node if one is returned
        if monopheltic_outgroup is not None:
            tree.set_outgroup(monopheltic_outgroup)
            rooted = True

        # returns the tree and if it was rooted
        return tree, rooted
    except Exception as e:
        print(f"Error rooting tree at outgroup: {e}")
        return tree, rooted