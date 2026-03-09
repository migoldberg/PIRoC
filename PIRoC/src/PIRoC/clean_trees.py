"""
clean_trees.py
Writes clean trees by removing contaminants and collapsing low support nodes.
"""

import os
import sys
import numpy as np
from datetime import datetime
from collections import Counter, defaultdict
from ete3 import Tree
from Bio import SeqIO

from .support import collapse_low_support_nodes

def fasta_in_tree_dir(tree_dir, og_id):
    """
    Helper function to check if the tree directory has the fasta file corresponding to the tree.
    """
    common_fasta_suffixes = [".fa", ".fasta", ".fna", ".faa"]
    for suffix in common_fasta_suffixes:
        fasta_path = os.path.join(tree_dir, f"{og_id}{suffix}")
        if os.path.exists(fasta_path):
            return fasta_path
    return False

def clean_fasta(og_id, contaminants, fasta_path, clean_dir):
    """
    Cleans the fasta by removing contaminants.
    """
    clean_fasta = (sequence for sequence in SeqIO.parse(fasta_path, "fasta") if sequence.id not in contaminants)
    SeqIO.write(clean_fasta, os.path.join(clean_dir, f"{og_id}.fa"), "fasta")

    return clean_fasta


def clean_trees(tree_dir, tree_suffix, sequence_classifications, output_dir, collapse_threshold, quiet):
    """
    Cleans trees by removing contaminants and collapsing low support nodes.
    """

    # creates a dictionary of trees
    trees = defaultdict(dict)

    # adds sequence ids and classifications to the trees dictionary
    for sequence_id, classification in sequence_classifications.items():
        og_id, sequence_name = sequence_id.split("::")
        trees[og_id][sequence_name] = classification

    # creates a directory for the clean orthogroups
    clean_dir = os.path.join(output_dir, "clean_orthogroups")
    os.makedirs(clean_dir, exist_ok=True)

    total_trees = len(trees)

    # loops through each tree in the trees dictionary
    for i, (og_id, sequences) in enumerate(trees.items(), 1):
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
        t.write(format=1, outfile=os.path.join(clean_dir, f"{og_id}{tree_suffix}"))

        if not quiet:
            if contaminants:
                print(f"[{og_id}] Removed {len(contaminants)} contaminant(s) and collapsed {nodes_collapsed} low support nodes (support < {collapse_threshold})")
            elif nodes_collapsed > 0:
                print(f"[{og_id}] Collapsed {nodes_collapsed} low support nodes (support < {collapse_threshold})")

            fasta_path = fasta_in_tree_dir(tree_dir, og_id)
            if fasta_path:
                clean_fasta(og_id, contaminants, fasta_path, clean_dir)
                print(f"[{og_id}] FASTA file cleaned")
            else:
                print(f'[{og_id}] No fasta file found in tree directory')
        else:
            fasta_path = fasta_in_tree_dir(tree_dir, og_id)
            if fasta_path:
                clean_fasta(og_id, contaminants, fasta_path, clean_dir)
            print(f"\r Cleaning Trees: [{i}/{total_trees}]", end="", flush=True)

    if quiet:
        print()

    return clean_dir