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
    trees = defaultdict(dict)
    for sequence_id, classification in sequence_classifications.items():
        og_id, sequence_name = sequence_id.split("::")
        trees[og_id][sequence_name] = classification

    clean_tree_dir = os.path.join(output_dir, "clean_trees")
    os.makedirs(clean_tree_dir, exist_ok=True)

    for og_id, sequences in trees.items():
        contaminants = [name for name, cls in sequences.items() if cls == "CONTAMINANT"]

        tree_path = os.path.join(tree_dir, f"{og_id}{tree_suffix}")
        t = Tree(tree_path, format=0)

        for contaminant in contaminants:
            matching_leaves = t.get_leaves_by_name(contaminant)
            for leaf in matching_leaves:
                leaf.detach()

        t, nodes_collapsed = collapse_low_support_nodes(t, collapse_threshold)
        t.write(format=1, outfile=os.path.join(clean_tree_dir, f"{og_id}{tree_suffix}"))

        if contaminants:
            print(f"[{og_id}] Removed {len(contaminants)} contaminant(s) and collapsed {nodes_collapsed} low support nodes")

        return len(contaminants)