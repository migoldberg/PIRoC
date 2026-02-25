"""
clean_trees.py
Writes clean trees by removing contaminants and collapsing low support nodes.
"""

import os
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree

from .support import collapse_low_support_nodes

def clean_trees(tree_dir, tree_suffix, sequence_classifications, output_dir, collapse_threshold):
    """
    Cleans trees by removing contaminants and collapsing low support nodes.
    """

    # creates a dictionary of trees
    trees = defaultdict(dict)

    # adds sequence ids and classifications to the trees dictionary
    for sequence_id, classification in sequence_classifications.items():
        og_id, sequence_name = sequence_id.split("::")
        trees[og_id][sequence_name] = classification

    # creates a directory for the clean trees
    clean_tree_dir = os.path.join(output_dir, "clean_trees")
    os.makedirs(clean_tree_dir, exist_ok=True)

    # loops throug heach tree in the trees dictionary
    for og_id, sequences in trees.items():
        # creates a list of contaminants for the tree
        contaminants = [name for name, cls in sequences.items() if cls == "CONTAMINANT"]

        # gets the path to the tree
        tree_path = os.path.join(tree_dir, f"{og_id}{tree_suffix}")

        # loads the tree using ete3
        t = Tree(tree_path, format=0)

        # loops through each contaminant in the tree
        for contaminant in contaminants:
            # finds the leaves in the tree that match the contaminant name
            matching_leaves = t.get_leaves_by_name(contaminant)

            # loops through each contaminant leaf in the tree and detaches it from the tree
            for leaf in matching_leaves:
                leaf.detach()

        # collapses low support nodes in the tree as done during the PIRoC run
        t, nodes_collapsed = collapse_low_support_nodes(t, collapse_threshold)

        # writes the clean tree to a new file
        t.write(format=1, outfile=os.path.join(clean_tree_dir, f"{og_id}{tree_suffix}"))

        # prints how many contaminants were removed and how many nodes were collapsed to the console
        if contaminants:
            print(f"[{og_id}] Removed {len(contaminants)} contaminant(s) and collapsed {nodes_collapsed} low support nodes")

        # returns the number of contaminants removed
        return len(contaminants)