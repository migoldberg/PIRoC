"""
support.py
Provides functions for finding the bootstrap support of a node, and to collapse low support nodes in the tree.
"""

import os 
import re
import sys
import argparse
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

def collapse_low_support_nodes(tree, threshold):
    if threshold is None:
        threshold = 50.0
    
    nodes_collapsed = 0

    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue

        support = getattr(node, 'support', None)

        if support is not None and 0.0 <= support <= 1.0:
            support = support * 100.0

        if support is not None and support < threshold:
            node.delete()
            nodes_collapsed += 1

    return nodes_collapsed

def get_node_bootstrap(node):
    if hasattr(node, "support") and node.support is not None:
        bs = float(node.support)
        if 0.0 <= bs <= 1.0:
            return bs * 100.0
        return bs
    else:
        return 0