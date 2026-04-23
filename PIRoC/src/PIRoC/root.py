"""
root.py
Attempts to root the tree at a contaminant clade if possible.
"""

import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree
from .metadata import parse_species_name, get_group_from_species_name

def find_monophyletic_contaminant(tree, contaminant_leaves):
    """
    Tries to find the most monophyletic clade of contaminant sequences in the tree.
    """

    # creates a set of the contaminant leaves
    contaminant_set = set(contaminant_leaves)

    # the most monophyletic contaminant node
    monophyletic_contaminant_node = None

    # the number of sequences in the most monophyletic contaminant node
    monophyletic_contaminant_node_leaves = 0

    # loops throug each node in the tree in postorder
    for idx, node in enumerate(tree.traverse("postorder")):
        # adds an index feature to the node
        node.add_feature("index", idx)

        if node.is_leaf():
            # if the node is a leaf then there is one leaf
            num_leaves = 1

            # if the node is a contaminant leaf then there is one contaminant leaf
            num_contaminant_leaves = 1 if node in contaminant_set else 0
        else:
            # if the node is not a leaf then there are multiple leaves
            num_leaves = sum(child.num_leaves for child in node.children)

            # finds the number of contaminant leaves which are children of the node
            num_contaminant_leaves = sum(child.num_contaminant_leaves for child in node.children)

        if node.up is not None:
            # if the node has a parent then add the number of leaves and contaminant leaves to the node features
            node.add_feature("num_leaves", num_leaves)
            node.add_feature("num_contaminant_leaves", num_contaminant_leaves)

            # calculates the fraction of leaves in the node that are contaminants
            contaminant_fraction = (num_contaminant_leaves / num_leaves) if num_leaves else 0.0

            # adds the contaminant fraction to the node features
            node.add_feature("contaminant_fraction", contaminant_fraction)

            # checks if the node is monophyletic
            is_monophyletic = (num_contaminant_leaves == num_leaves) and (num_contaminant_leaves > 0)

            # if the node is more monophyletic than the current most monophyletic contaminant node then update the most monophyletic contaminant node
            if is_monophyletic and num_contaminant_leaves > monophyletic_contaminant_node_leaves:
                monophyletic_contaminant_node_leaves = num_contaminant_leaves
                monophyletic_contaminant_node = node

    # returns the most monophyletic contaminant node
    return monophyletic_contaminant_node

def root_tree_at_contaminant(tree, contaminant_leaves):
    """
    Attempts to root the tree at a monophyletic contaminant clade.
    """

    # assumes the tree is not rooted (by this script, if this does not re-root it than the tree-builder may have rooted it somewhere)
    rooted = False

    # if there are no contaminant leaves then return the tree
    if not contaminant_leaves:
        return tree, rooted

    # if there is only one contaminant leaf then root the tree there and return the tree
    if len(contaminant_leaves) == 1:
        contaminant_leaf = contaminant_leaves[0]
        tree.set_outgroup(contaminant_leaf)
        rooted = True
        return tree, rooted

    try:
        # finds the most monophyletic contaminant clade
        monophyletic_contaminant = find_monophyletic_contaminant(tree, contaminant_leaves)

        # root the tree at the most monophyletic contaminant node if one is returned (ete3 API: set_outgroup)
        if monophyletic_contaminant is not None:
            tree.set_outgroup(monophyletic_contaminant)
            rooted = True

        # returns the tree and if it was rooted
        return tree, rooted
    except Exception as e:
        print(f"Error rooting tree at contaminant: {e}")
        return tree, rooted