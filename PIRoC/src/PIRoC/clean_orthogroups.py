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
from .metadata import get_group_from_species_name

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

def clean_fasta(og_id, contaminants, fasta_path, clean_dir, species_to_group, contaminant_names):
    """
    Cleans the fasta by removing contaminants.
    """
    clean_fasta = (
        sequence for sequence in SeqIO.parse(fasta_path, "fasta")
        if sequence.id not in contaminants
        and get_group_from_species_name(sequence.id, species_to_group) not in contaminant_names
    )
    SeqIO.write(clean_fasta, os.path.join(clean_dir, f"{og_id}.fa"), "fasta")

    return clean_fasta

def clean_orthogroups(tree_dir, tree_suffix, sequence_classifications, contaminant_names, species_to_group, output_dir, loud):
    trees = defaultdict(dict)

    for sequence_id, classification in sequence_classifications.items():
        og_id, sequence_name = sequence_id.split("::")
        trees[og_id][sequence_name] = classification

    clean_dir = os.path.join(output_dir, "clean_orthogroups")
    os.makedirs(clean_dir, exist_ok=True)

    total_trees = len(trees)

    for i, (og_id, sequences) in enumerate(trees.items(), 1):
        contaminants = [name for name, cls in sequences.items() if cls == "CONTAMINANT"]
      
        fasta_path = fasta_in_tree_dir(tree_dir, og_id)

        if loud:
            if fasta_path:
                clean_fasta(og_id, contaminants, fasta_path, clean_dir, species_to_group, contaminant_names)
                print(f"[{og_id}] FASTA file cleaned")
            else:
                print(f'[{og_id}] No fasta file found in tree directory')
        else:
            fasta_path = fasta_in_tree_dir(tree_dir, og_id)
            if fasta_path:
                clean_fasta(og_id, contaminants, fasta_path, clean_dir, species_to_group, contaminant_names)
            print(f"\r  Cleaning Orthogroups: [{i}/{total_trees}]", end="", flush=True)
       
    return clean_dir
       