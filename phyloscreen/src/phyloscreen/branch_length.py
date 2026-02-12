"""
branch_length.py
Provides functions for identifying sequences occuring on long branches in the tree.
"""

import os 
import re
import sys
import argparse
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

def compute_branch_length_stats(tree):
    branch_lengths = []
    root = tree.get_tree_root()

    if not root:
        return None

    for leaf in tree.iter_leaves():
        branch_lengths.append(root.get_distance(leaf))

    if not branch_lengths:
        return None

    branch_distances = np.array(branch_lengths)
    mean = branch_distances.mean()
    std = branch_distances.std()

    long_branch_threshold = mean + (2 * std)

    return {
        "mean": mean,
        "std": std,
        "long_branch_threshold": long_branch_threshold
    }

def is_on_long_branch(leaf, tree, branch_length_stats):
    is_long_branch = False

    if branch_length_stats is None:
        return is_long_branch

    threshold = branch_length_stats['long_branch_threshold']
    mean_branch_length = branch_length_stats['mean']

    root = tree.get_tree_root()
    if not root:
        return is_long_branch
    
    leaf_root_to_tip_distance = root.get_distance(leaf)

    if leaf_root_to_tip_distance > threshold:
        is_long_branch = True
    else:
        is_long_branch = False

    return is_long_branch

    