"""
support.py
Provides functions for finding the bootstrap support of a node, and to collapse low support nodes in the tree.
"""

import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

def collapse_low_support_nodes(tree, threshold):
    """
    ** No longer used because trees are not outputted anymore so collapsing would not be helpful.

    Collapses low support nodes in the tree.
    """

    # if no threshold is provided then the tree is returned unchanged
    if threshold is None:
        return tree, 0

    # initializes a list of low support nodes
    low_support_nodes = []

    # loops through each node in the tree in postorder
    for node in tree.traverse("postorder"):
        # skips leaves and the root node
        if node.is_leaf() or node.is_root():
            continue
    
        # gets the bootstrap support of the node
        support = get_node_bootstrap(node)

        # if the bootstrap support is less than the threshold, the node is added to the list of low support nodes
        if support < threshold:
            low_support_nodes.append(node)
    
    # loops through each low support node and deletes it
    for node in low_support_nodes:
        node.delete()

    # returns the tree and the number of low support nodes collapsed
    return tree, len(low_support_nodes)

def get_node_bootstrap(node):
    """
    Gets the bootstrap support of a node.
    """

    # if the node has a support attribute and it is not None
    if hasattr(node, "support") and node.support is not None:
        # converts the support to a float
        bs = float(node.support)

        # if the bootstrap support is between 0 and 1 then normalize it by multiplying by 100
        if 0.0 <= bs <= 1.0:
            return bs * 100.0

        # return the bootstrap support
        return bs
    else:
        # if the node does not have a support attribute or it is None then return 0
        return 0