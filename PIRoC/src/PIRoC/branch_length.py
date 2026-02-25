"""
branch_length.py
Provides functions for identifying sequences occuring on long branches in the tree.
"""

import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

def compute_branch_length_stats(tree):
    """
    Computes the root-to-tip distance for each leaf in the tree and calcualtes their mean.
    Two standard deviations above the mean is used as a threshold for identifying sequences occuring on long branches.
    """
    # initializes an empty list to store the branch lengths
    branch_lengths = []
    # gets the root of the tree
    root = tree.get_tree_root()

    # if the tree is not rooted, None is returned
    if not root:
        return None

    # loops through each leaf in the tree and appends the root-to-tip distance to the list
    for leaf in tree.iter_leaves():
        branch_lengths.append(root.get_distance(leaf))

    if not branch_lengths:
        return None

    # converts the list of branch lengths to a numpy array
    branch_distances = np.array(branch_lengths)

    # calculates the mean of the branch lengths
    mean = branch_distances.mean()

    # calculates the standard deviation of the branch lengths
    std = branch_distances.std()

    # calculates the threshold for identifying sequences occuring on long branches (currently 2 std deviations above the mean)
    long_branch_threshold = mean + (2 * std)

    # returns the branch length statistics
    return {
        "mean": mean,
        "std": std,
        "long_branch_threshold": long_branch_threshold
    }

def is_on_long_branch(leaf, tree, branch_length_stats):
    """
    Determines if a sequence is occuring on a long branch in the tree.
    """
    # starts assuming the sequence is not on a long branch
    is_long_branch = False

    # if the branch length statistics do not exist, returns False
    if branch_length_stats is None:
        return is_long_branch

    # gets the long branch threshold from the branch length stats
    threshold = branch_length_stats['long_branch_threshold']

    # gets the root of the tree
    root = tree.get_tree_root()

    # if the tree is not rooted, returns False since the root-to-tip distance cannot be calculated
    if not root:
        return is_long_branch
    
    # calculates the root-to-tip distance for the sequence
    leaf_root_to_tip_distance = root.get_distance(leaf)

    # if the root-to-tip distance is greater than the threshold, the sequence is on a long branch
    if leaf_root_to_tip_distance > threshold:
        is_long_branch = True
    else:
        is_long_branch = False

    return is_long_branch

    