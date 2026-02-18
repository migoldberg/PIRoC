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
    outgroup_set = set(outgroup_leaves)

    monophyletic_outgroup_node = None
    monophletic_outgroup_node_leaves = 0

    for idx, node in enumerate(tree.traverse("postorder")):
        node.add_feature("index", idx)
        if node.is_leaf():
            num_leaves = 1
            num_outgroup_leaves = 1 if node in outgroup_set else 0
        else:
            num_leaves = sum(child.num_leaves for child in node.children)
            num_outgroup_leaves = sum(child.num_outgroup_leaves for child in node.children)
        if node.up is not None:
            node.add_feature("num_leaves", num_leaves)
            node.add_feature("num_outgroup_leaves", num_outgroup_leaves)

            outgroup_fraction = (num_outgroup_leaves / num_leaves) if num_leaves else 0.0
            node.add_feature("outgroup_fraction", outgroup_fraction)

            is_monophyletic = (num_outgroup_leaves == num_leaves) and (num_outgroup_leaves > 0)

            if is_monophyletic and num_outgroup_leaves > monophletic_outgroup_node_leaves:
                monophletic_outgroup_node_leaves = num_outgroup_leaves
                monophyletic_outgroup_node = node

    return monophyletic_outgroup_node


def root_tree_at_outgroup(tree, outgroup_leaves):
    rooted = False

    if not outgroup_leaves:
        return tree, rooted

    if len(set(outgroup_leaves)) < len(outgroup_leaves):
        return tree, rooted

    if len(outgroup_leaves) == 1:
        outgroup_leaf = outgroup_leaves[0]
        tree.set_outgroup(outgroup_leaf)
        rooted = True
        return tree, rooted

    try:
        monopheltic_outgroup = find_monophyletic_outgroup(tree, outgroup_leaves)

        if monopheltic_outgroup is not None:
            tree.set_outgroup(monopheltic_outgroup)
            rooted = True

        return tree, rooted
    except Exception as e:
        print(f"Error rooting tree at outgroup: {e}")
        return tree, rooted