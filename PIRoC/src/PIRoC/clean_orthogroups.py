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

from .metadata import get_group_from_species_name, parse_species_name

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

def write_removed_fasta(contaminants, fasta_path, removed_dir):
    """
    Writes removed contaminant sequences to an overall fasta and to one fasta
    per species inside removed_dir. Returns the overall fasta path.
    """
    overall_fasta_path = os.path.abspath(os.path.join(removed_dir, "all_contaminants.fa"))
    if not contaminants:
        return overall_fasta_path

    contaminant_ids = set(contaminants)
    removed_sequences = [
        sequence for sequence in SeqIO.parse(fasta_path, "fasta")
        if sequence.id in contaminant_ids
    ]
    if not removed_sequences:
        return overall_fasta_path

    os.makedirs(removed_dir, exist_ok=True)
    species_dir = os.path.join(removed_dir, "by_species")
    os.makedirs(species_dir, exist_ok=True)

    with open(overall_fasta_path, "a") as handle:
        SeqIO.write(removed_sequences, handle, "fasta")

    sequences_by_species = defaultdict(list)
    for sequence in removed_sequences:
        species = parse_species_name(sequence.id) or "UNKNOWN"
        sequences_by_species[species].append(sequence)

    for species, sequences in sequences_by_species.items():
        species_fasta_path = os.path.join(species_dir, f"{species}.fa")
        with open(species_fasta_path, "a") as handle:
            SeqIO.write(sequences, handle, "fasta")

    return overall_fasta_path

def clean_fasta(og_id, contaminants, fasta_path, clean_dir, species_to_group, contaminant_names, keep_known_contaminants):
    """
    Cleans the fasta by removing contaminants.
    """
    clean_fasta = (
        sequence for sequence in SeqIO.parse(fasta_path, "fasta")
        if sequence.id not in contaminants
        and (
            keep_known_contaminants
            or get_group_from_species_name(sequence.id, species_to_group) not in contaminant_names
        )
    )
    SeqIO.write(clean_fasta, os.path.join(clean_dir, f"{og_id}.fa"), "fasta")

    return clean_fasta

def clean_orthogroups(tree_dir, tree_suffix, sequence_classifications, contaminant_names, species_to_group, output_dir, loud, keep_known_contaminants, remove_flags=False):
    trees = defaultdict(dict)

    for sequence_id, classification in sequence_classifications.items():
        og_id, sequence_name = sequence_id.split("::")
        trees[og_id][sequence_name] = classification

    clean_dir = os.path.join(output_dir, "clean_orthogroups")
    os.makedirs(clean_dir, exist_ok=True)

    removed_dir = os.path.join(output_dir, "removed_contaminants")
    species_dir = os.path.join(removed_dir, "by_species")
    os.makedirs(species_dir, exist_ok=True)

    overall_fasta_path = os.path.join(removed_dir, "all_contaminants.fa")
    with open(overall_fasta_path, "w"):
        pass
    for existing in os.listdir(species_dir):
        existing_path = os.path.join(species_dir, existing)
        if os.path.isfile(existing_path) and existing.endswith(".fa"):
            os.remove(existing_path)

    total_trees = len(trees)

    removal_classes = {"CONTAMINANT"}
    if remove_flags:
        removal_classes.add("FLAG")

    for i, (og_id, sequences) in enumerate(trees.items(), 1):
        contaminants = [name for name, cls in sequences.items() if cls in removal_classes]
      
        fasta_path = fasta_in_tree_dir(tree_dir, og_id)

        if fasta_path:
            write_removed_fasta(contaminants, fasta_path, removed_dir)
            clean_fasta(og_id, contaminants, fasta_path, clean_dir, species_to_group, contaminant_names, keep_known_contaminants)
            if loud:
                print(f"[{og_id}] FASTA file cleaned")
        elif loud:
            print(f'[{og_id}] No fasta file found in tree directory')

        if not loud:
            print(f"\r  Cleaning Orthogroups: [{i}/{total_trees}]", end="", flush=True)

    return clean_dir
       