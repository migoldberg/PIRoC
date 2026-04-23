"""
support.py
Provides functions for finding the bootstrap support of a node.
"""

import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

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