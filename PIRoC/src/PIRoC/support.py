"""
support.py
Provides functions for finding the bootstrap support of a node, and to collapse low support nodes in the tree.
"""

import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

def collapse_low_support_nodes(tree, threshold):
    if threshold is None:
        return tree, 0

    low_support_nodes = []

    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
    
        support = get_node_bootstrap(node)

        if support < threshold:
            low_support_nodes.append(node)
    
    for node in low_support_nodes:
        node.delete()

    return tree, len(low_support_nodes)

def get_node_bootstrap(node):
    if hasattr(node, "support") and node.support is not None:
        bs = float(node.support)
        if 0.0 <= bs <= 1.0:
            return bs * 100.0
        return bs
    else:
        return 0